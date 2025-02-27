/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2014, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
 * THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the Open Source Modelica
 * Consortium (OSMC) Public License (OSMC-PL) are obtained
 * from OSMC, either from the above address,
 * from the URLs: http://www.ida.liu.se/projects/OpenModelica or
 * http://www.openmodelica.org, and in the OpenModelica distribution.
 * GNU version 3 is obtained from: http://www.gnu.org/copyleft/gpl.html.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of  MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

encapsulated package ExpressionSolve
" file:        ExpressionSolve.mo
  package:     ExpressionSolve
  description: ExpressionSolve


  This file contains the module ExpressionSolve, which contains functions
  to solve a DAE.Exp for a DAE.Exp"

// public imports
public import Absyn;
public import AbsynUtil;
public import DAE;

// protected imports
protected import ComponentReference;
protected import Debug;
protected import Differentiate;
protected import ElementSource;
protected import Expression;
protected import ExpressionDump;
protected import ExpressionSimplify;
protected import Flags;
protected import List;
protected import Inline;
protected import BackendDAEUtil;
protected import BackendEquation;
protected import BackendVariable;

// =============================================================================
// section for postOptModule >>solveSimpleEquations<<
//
// solve simple equations otherwise detect EQUATIONSYSTEM
// =============================================================================

public function solveSimpleEquations
  input output BackendDAE.BackendDAE dae;
algorithm
  dae.eqs := list(
    match syst
      local
        BackendDAE.StrongComponents comps;
        array<Integer> ass1 "eqn := ass1[var]";
        array<Integer> ass2 "var := ass2[eqn]";

      case BackendDAE.EQSYSTEM(matching = BackendDAE.MATCHING(comps=comps, ass1=ass1, ass2=ass2))
        algorithm
          comps := list(
            match comp
              local
                BackendDAE.Equation eqn;
                BackendDAE.Var var;
                Integer eindex, vindx;
                Boolean solved;
                BackendDAE.StrongComponent tmpComp;

              case BackendDAE.SINGLEEQUATION() algorithm
                BackendDAE.SINGLEEQUATION(eqn=eindex, var=vindx) := comp;
                eqn := BackendEquation.get(syst.orderedEqs, eindex);
                tmpComp := comp;
                if BackendEquation.isEquation(eqn) then
                  var := BackendVariable.getVarAt(syst.orderedVars, vindx);
                  (eqn, solved) := solveSimpleEquation(eqn, var, dae.shared);
                  syst.orderedEqs := BackendEquation.setAtIndex(syst.orderedEqs, eindex, eqn);
                  if not solved then
                    tmpComp := BackendDAE.EQUATIONSYSTEM({eindex}, {vindx}, BackendDAE.EMPTY_JACOBIAN(), BackendDAE.JAC_NONLINEAR(), false);
                  end if;
                end if;
              then tmpComp;

              else comp;
            end match
          for comp in comps);
          syst.matching := BackendDAE.MATCHING(ass1, ass2, comps);
        then syst;

      else syst;
    end match
  for syst in dae.eqs);
end solveSimpleEquations;

protected function solveSimpleEquation
  input output BackendDAE.Equation eqn;
  input BackendDAE.Var var "solve eqn with respect to var";
  input BackendDAE.Shared shared;
  output Boolean solved;
protected
  DAE.ComponentRef cr;
  DAE.Exp e1,e2,varexp,e;
  BackendDAE.EquationAttributes attr;
  DAE.ElementSource source;
  Boolean isContinuousIntegration = BackendDAEUtil.isSimulationDAE(shared);
algorithm
  BackendDAE.EQUATION(exp=e1, scalar=e2, source=source, attr=attr) := eqn;
  BackendDAE.VAR(varName = cr) := var;
  varexp := Expression.crefExp(cr);
  if BackendVariable.isStateVar(var) then
    varexp := Expression.expDer(varexp);
    cr := ComponentReference.crefPrefixDer(cr);
  end if;

  // phi: aren't the types of e1 and e2 the same? Can we make the equation have a type?
  if (Types.isIntegerOrRealOrSubTypeOfEither(Expression.typeof(e1)) and Types.isIntegerOrRealOrSubTypeOfEither(Expression.typeof(e2))) then
    (e1, e2) := preprocessingSolve(e1, e2, varexp, NONE(), SOME(shared.functionTree), NONE(), 0,  false);
  end if;

  try
    e := solve2(e1, e2, varexp, SOME(shared.functionTree), NONE(), false, isContinuousIntegration);
    source := ElementSource.addSymbolicTransformationSolve(true, source, cr, e1, e2, e, {});
    eqn := BackendEquation.generateEquation(varexp, e, source, attr);
    solved := true;
  else
    // only return new eqn if it can be solved explicitely because intermediate results can be numerically bad
    // solves ticket #4293
    // ToDo: do other preprocessing like multiplying by divisors?
    solved := false;
  end try;
end solveSimpleEquation;

protected function printTryToSolve
  "for debugging"
  input String instanceName "getInstanceName from caller";
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
algorithm
  print(instanceName + " tries to solve: " +
    ExpressionDump.printExpStr(inExp1) + " = " + ExpressionDump.printExpStr(inExp2) +
    "\nwith respect to: " + ExpressionDump.printExpStr(inExp3) + "\n");
end printTryToSolve;

public function solve
"Solves an equation consisting of a right hand side (rhs) and a
  left hand side (lhs), with respect to the expression given as
  third argument, usually a variable."
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  input Option<DAE.FunctionTree> functions = NONE() "need for solve modelica functions";
  output DAE.Exp outExp;
  output list<DAE.Statement> outAsserts;
protected
  list<BackendDAE.Equation> dummy1;
  list<DAE.ComponentRef> dummy2;
  Integer dummyI;
algorithm
  //printTryToSolve(getInstanceName(), inExp1, inExp2, inExp3);

  (outExp,outAsserts,dummy1, dummy2, dummyI) := matchcontinue inExp1
    case _ then solveSimple(inExp1, inExp2, inExp3, 0);
    case _ then solveSimple(inExp2, inExp1, inExp3, 0);
    case _ then solveWork(inExp1, inExp2, inExp3, NONE(), functions, NONE(), 0, false, false);
    else equation
      if Flags.isSet(Flags.FAILTRACE) then
        Error.addInternalError("Failed to solve \"" + ExpressionDump.printExpStr(inExp1) + " = " + ExpressionDump.printExpStr(inExp2) + "\" w.r.t. \"" + ExpressionDump.printExpStr(inExp3) + "\"", sourceInfo());
      end if;
    then fail();
  end matchcontinue;

 (outExp,_) := ExpressionSimplify.simplify1(outExp);
end solve;


public function solve2
"Solves an equation with modelica function consisting of a right hand side (rhs) and a
  left hand side (lhs), with respect to the expression given as
  third argument, usually a variable.
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  input Option<DAE.FunctionTree> functions "need for solve modelica functions";
  input Option<Integer> uniqueEqIndex "offset for tmp vars";
  input Boolean doInline = true;
  input Boolean isContinuousIntegration = false;
  output DAE.Exp outExp;
  output list<DAE.Statement> outAsserts;
  output list<BackendDAE.Equation> eqnForNewVars "eqn for tmp vars";
  output list<DAE.ComponentRef> newVarsCrefs;
protected
  Integer dummyI;
algorithm
  //printTryToSolve(getInstanceName(), inExp1, inExp2, inExp3);

  (outExp,outAsserts,eqnForNewVars,newVarsCrefs,dummyI) := matchcontinue inExp1
    case _ then solveSimple(inExp1, inExp2, inExp3, 0);
    case _ then solveSimple(inExp2, inExp1, inExp3, 0);
    case _ then solveWork(inExp1, inExp2, inExp3, NONE(), functions, uniqueEqIndex, 0, doInline, isContinuousIntegration);
    else equation
      if Flags.isSet(Flags.FAILTRACE) then
        Error.addInternalError("Failed to solve \"" + ExpressionDump.printExpStr(inExp1) + " = " + ExpressionDump.printExpStr(inExp2) + "\" w.r.t. \"" + ExpressionDump.printExpStr(inExp3) + "\"", sourceInfo());
      end if;
    then fail();
  end matchcontinue;
end solve2;


protected function solveWork
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  input Option<DAE.Exp> optCond "condition from an if expression";
  input Option<DAE.FunctionTree> functions;
  input Option<Integer> uniqueEqIndex "offset for tmp vars";
  input Integer idepth;
  input Boolean doInline;
  input Boolean isContinuousIntegration;
  output DAE.Exp outExp;
  output list<DAE.Statement> outAsserts;
  output list<BackendDAE.Equation> eqnForNewVars "eqn for tmp vars";
  output list<DAE.ComponentRef> newVarsCrefs;
  output Integer depth;
protected
  DAE.Exp e1, e2;
  list<BackendDAE.Equation> eqnForNewVars1, eqnForNewVars2;
  list<DAE.ComponentRef> newVarsCrefs1, newVarsCrefs2;
algorithm
  (e1, e2, eqnForNewVars1, newVarsCrefs1, depth) := matchcontinue inExp1
    case _ then preprocessingSolve(inExp1, inExp2, inExp3, optCond, functions, uniqueEqIndex, idepth, doInline);
    else
      equation
        if Flags.isSet(Flags.FAILTRACE) then
          Debug.trace("\n-ExpressionSolve.preprocessingSolve failed:\n");
          Debug.trace(ExpressionDump.printExpStr(inExp1) + " = " + ExpressionDump.printExpStr(inExp2));
          Debug.trace(" with respect to: " + ExpressionDump.printExpStr(inExp3));
        end if;
      then (inExp1,inExp2,{},{}, idepth);
  end matchcontinue;

  (outExp, outAsserts, eqnForNewVars2, newVarsCrefs2, depth) := matchcontinue e1
    case _ then solveIfExp(e1, e2, inExp3, optCond, functions, uniqueEqIndex, depth, doInline, isContinuousIntegration);
    case _ then solveSimple(e1, e2, inExp3, depth);
    case _ then solveLinearSystem(e1, e2, inExp3, functions, depth);
   end matchcontinue;

  eqnForNewVars := listAppend(eqnForNewVars1, eqnForNewVars2);
  newVarsCrefs := listAppend(newVarsCrefs1, newVarsCrefs2);
end solveWork;

protected function solveSimple
"Solves simple equations like
  a = f(..)
  der(a) = f(..)
  -a = f(..)
  -der(a) = f(..)"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  input Integer idepth;
  output DAE.Exp outExp;
  output list<DAE.Statement> outAsserts;
  output list<BackendDAE.Equation> eqnForNewVars = {} "eqn for tmp vars";
  output list<DAE.ComponentRef> newVarsCrefs = {};
  output Integer odepth = idepth;

algorithm
  //printTryToSolve(getInstanceName(), inExp1, inExp2, inExp3);

  (outExp,outAsserts) := match (inExp1,inExp2,inExp3)
    local
      DAE.ComponentRef cr,cr1;
      DAE.Type tp;
      DAE.Exp e1,e2,res,e11;
      Real r, r2;
      list<DAE.Statement> asserts;

    // special case when already solved, cr1 = rhs, otherwise division by zero when dividing with derivative
    case (DAE.CREF(componentRef = cr1),_,DAE.CREF(componentRef = cr))
      guard ComponentReference.crefEqual(cr, cr1) and (not Expression.expHasCrefNoPreOrStart(inExp2, cr))
      then
        (inExp2,{});
    case (DAE.CALL(path = Absyn.IDENT(name = "der"),expLst = {DAE.CREF(componentRef = cr1)}),_,DAE.CALL(path = Absyn.IDENT(name = "der"),expLst = {DAE.CREF(componentRef = cr)}))
      guard  ComponentReference.crefEqual(cr, cr1) and (not Expression.expHasDerCref(inExp2, cr))
      then
        (inExp2,{});

    // -cr = exp
    case (DAE.UNARY(operator = DAE.UMINUS(), exp = DAE.CREF(componentRef = cr1)),_,DAE.CREF(componentRef = cr))
      guard ComponentReference.crefEqual(cr1,cr) and (not Expression.expHasCrefNoPreOrStart(inExp2,cr))
      then
        (Expression.negate(inExp2),{});
    case (DAE.UNARY(operator = DAE.UMINUS_ARR(), exp = DAE.CREF(componentRef = cr1)),_,DAE.CREF(componentRef = cr))
      guard ComponentReference.crefEqual(cr1,cr) and (not Expression.expHasCrefNoPreOrStart(inExp2,cr)) // cr not in e2
      then
        (Expression.negate(inExp2),{});
    case (DAE.UNARY(operator = DAE.UMINUS(), exp = DAE.CALL(path = Absyn.IDENT(name = "der"),expLst = {DAE.CREF(componentRef = cr1)})),_,DAE.CALL(path = Absyn.IDENT(name = "der"),expLst = {DAE.CREF(componentRef = cr)}))
      guard ComponentReference.crefEqual(cr1,cr) and (not  Expression.expHasDerCref(inExp2,cr)) // cr not in e2
      then
        (Expression.negate(inExp2),{});
    case (DAE.UNARY(operator = DAE.UMINUS_ARR(), exp = DAE.CALL(path = Absyn.IDENT(name = "der"),expLst = {DAE.CREF(componentRef = cr1)})),_,DAE.CALL(path = Absyn.IDENT(name = "der"),expLst = {DAE.CREF(componentRef = cr)}))
      guard ComponentReference.crefEqual(cr1,cr) and (not Expression.expHasDerCref(inExp2,cr))
      then
        (Expression.negate(inExp2),{});

    // !cr = exp
    case (DAE.LUNARY(operator = DAE.NOT(), exp = DAE.CREF(componentRef = cr1)),_,DAE.CREF(componentRef = cr))
      guard ComponentReference.crefEqual(cr1,cr) and (not Expression.expHasCrefNoPreOrStart(inExp2,cr))
      then
        (Expression.negate(inExp2),{});

    // Integer(enumcr) = ...
    case (DAE.CALL(path = Absyn.IDENT(name = "Integer"),expLst={DAE.CREF(componentRef = cr1)}),_,DAE.CREF(componentRef = cr,ty=tp))
      guard ComponentReference.crefEqual(cr, cr1) and (not  Expression.expHasCrefNoPreorDer(inExp2,cr))
      equation
        asserts = generateAssertType(tp,cr,inExp3,{});
      then (DAE.CAST(tp,inExp2),asserts);
      else fail();
  end match;
end solveSimple;

protected function generateAssertType
  input DAE.Type tp;
  input DAE.ComponentRef cr;
  input DAE.Exp iExp;
  input list<DAE.Statement> inAsserts;
  output list<DAE.Statement> outAsserts;
algorithm
  outAsserts := match(tp,cr,iExp,inAsserts)
    local
      Absyn.Path path,p1,pn;
      list<String> names;
      Integer n;
      DAE.Exp e1,en,e,es;
      String s1,sn,se,estr,crstr;
    case (DAE.T_ENUMERATION(path=path,names=names),_,_,_)
      equation
        p1 = AbsynUtil.suffixPath(path,listHead(names));
        e1 = DAE.ENUM_LITERAL(p1,1);
        n = listLength(names);
        pn = AbsynUtil.suffixPath(path,listGet(names,n));
        en = DAE.ENUM_LITERAL(p1,n);
        s1 = AbsynUtil.pathString(p1);
        sn = AbsynUtil.pathString(pn);
        _ = ExpressionDump.printExpStr(iExp);
        crstr = ComponentReference.printComponentRefStr(cr);
        estr = "Expression for " + crstr + " out of min(" + s1 + ")/max(" + sn + ") = ";
        // iExp >= e1 and iExp <= en
        e = DAE.LBINARY(DAE.RELATION(iExp,DAE.GREATEREQ(DAE.T_INTEGER_DEFAULT),e1,-1,NONE()),DAE.AND(DAE.T_BOOL_DEFAULT),
                                     DAE.RELATION(iExp,DAE.LESSEQ(DAE.T_INTEGER_DEFAULT),en,-1,NONE()));
        es = Expression.makePureBuiltinCall("String", {iExp,DAE.SCONST("d")}, DAE.T_STRING_DEFAULT);
        es = DAE.BINARY(DAE.SCONST(estr),DAE.ADD(DAE.T_STRING_DEFAULT),es);
      then
        DAE.STMT_ASSERT(e,es,DAE.ASSERTIONLEVEL_ERROR,DAE.emptyElementSource)::inAsserts;
    else inAsserts;
  end match;
end generateAssertType;

public function preprocessingSolve
"
 preprocessing for solve1,
 sorting and split terms , with respect to the expression given as
 third argument.

 {f(x,y), g(x,y),x} -> {h(x), k(y)}

 author: Vitalij Ruge
"

  input output DAE.Exp x "lhs";
  input output DAE.Exp y "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  input Option<DAE.Exp> optCond "condition from an if expression";
  input Option<DAE.FunctionTree> functions;
  input Option<Integer> uniqueEqIndex "offset for tmp vars";
  input Integer idepth;
  input Boolean doInline;
  output list<BackendDAE.Equation> eqnForNewVars = {} "eqn for tmp vars";
  output list<DAE.ComponentRef> newVarsCrefs = {};
  output Integer depth = idepth;

 protected
  DAE.Exp res;
  list<DAE.Exp> lhs, rhs;
  list<DAE.Exp> lhsWithX, rhsWithX, lhsWithoutX, rhsWithoutX, eWithX, factorWithX, factorWithoutX;
  DAE.Exp lhsX, rhsX, lhsY, rhsY, N;
  DAE.ComponentRef cr;
  DAE.Boolean con, new_x, inlineFun = true;
  Integer iter;
  Integer numSimplifed = 0 ;

 algorithm

   // split and sort
   (lhsX, lhsY) := preprocessingSolve5(x, inExp3,true);
   (rhsX, rhsY) := preprocessingSolve5(y, inExp3,true);
   x := Expression.expSub(lhsX, rhsX);
   y := Expression.expSub(rhsY, lhsY);

   con := not Expression.isCref(x);
   iter := 0;
   if con then
     x := unifyFunCalls(x, inExp3);
   end if;
   while con and iter < 1000 and (not Expression.isCref(x)) loop
     (x, y, con) := preprocessingSolve2(x,y, inExp3);
     (x, y, new_x) := preprocessingSolve3(x,y, inExp3);
     con := con or new_x;
     while new_x loop
       (x, y, new_x) := preprocessingSolve3(x,y, inExp3);
     end while;

     if Expression.isCref(x) then
       break;
     end if;
     (x, y, new_x) := removeSimpleCalls(x,y, inExp3);
     con := con or new_x;
     (x, y, new_x) := preprocessingSolve4(x,y, inExp3);
     con := new_x or con;
     // TODO: use new defined function, which missing in the cpp runtime
     if isSome(uniqueEqIndex) and not stringEqual(Config.simCodeTarget(), "Cpp") then
       (x, y, new_x, eqnForNewVars, newVarsCrefs, depth) := preprocessingSolveTmpVars(x, y, inExp3, optCond, Util.getOption(uniqueEqIndex), eqnForNewVars, newVarsCrefs, depth);
       con := new_x or con;
     end if;

     if (not con) then
       if (numSimplifed < 3) then
         (x, con) := ExpressionSimplify.simplify(x);
         numSimplifed := numSimplifed + 1;
       end if;
       // Z/N = rhs -> Z = rhs*N
       (x,N) := Expression.makeFraction(x);
       if not Expression.isOne(N) then
         //print("\nx ");print(ExpressionDump.printExpStr(x));print("\nN ");print(ExpressionDump.printExpStr(N));
         new_x := true;
         y := Expression.expMul(y,N);
       end if;

       con := new_x or con;
     end if;

     if con  then
       (lhsX, lhsY) := preprocessingSolve5(x, inExp3, true);
       (rhsX, rhsY) := preprocessingSolve5(y, inExp3, false);
       x := Expression.expSub(lhsX, rhsX);
       y := Expression.expSub(rhsY, lhsY);
     elseif doInline and inlineFun then
       iter := iter + 50;
       if inlineFun then
         (x,con) := solveFunCalls(x, inExp3, functions);
         inlineFun := false;
         if con then
           numSimplifed := 0;
         end if;
       end if;
     end if;

     iter := iter + 1;
     //print("\nx ");print(ExpressionDump.printExpStr(x));print("\ny ");print(ExpressionDump.printExpStr(y));
   end while;

   y := ExpressionSimplify.simplify1(y);

/*
   if not Expression.expEqual(inExp1,x) then
     print("\nIn: ");print(ExpressionDump.printExpStr(inExp1));print(" = ");print(ExpressionDump.printExpStr(inExp2));
     print("\nOut: ");print(ExpressionDump.printExpStr(x));print(" = ");print(ExpressionDump.printExpStr(y));
     print("\t w.r.t ");print(ExpressionDump.printExpStr(inExp3));
   end if;
*/
end preprocessingSolve;

protected function preprocessingSolve2
"
 helprer function for preprocessingSolve
 e.g.
   x/(x+c1) = -c2 --> x + (x+c1)*c2 = 0

 author: Vitalij Ruge
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";

  output DAE.Exp olhs;
  output DAE.Exp orhs;
  output Boolean con "continue";

algorithm

  (olhs, orhs, con) := match (inExp1)
    local
     DAE.Exp e,a, b, fb, fa, ga, lhs;
     DAE.Type tp;
     DAE.Operator op;
     list<DAE.Exp> eWithX, factorWithX, factorWithoutX;
     DAE.Exp pWithX, pWithoutX;

    // -f(a) = b => f(a) = -b
    case DAE.UNARY(op as DAE.UMINUS(), fa)
      guard expHasCref(fa, inExp3) and not expHasCref(inExp2, inExp3)
      equation
        b = Expression.negate(inExp2);
      then (fa, b, true);

    case DAE.UNARY(op as DAE.UMINUS_ARR(), fa)
      guard expHasCref(fa, inExp3) and not expHasCref(inExp2, inExp3)
      equation
        b = Expression.negate(inExp2);
      then (fa, b, true);

    // b/f(a) = rhs  => f(a) = b/rhs solve for a
    case DAE.BINARY(b,DAE.DIV(_),fa)
      guard expHasCref(fa, inExp3) and (not expHasCref(b, inExp3)) and (not expHasCref(inExp2, inExp3))
      equation
        e = Expression.makeDiv(b, inExp2);
      then(fa, e, true);

    // b*f(a) = rhs  => f(a) = rhs/b solve for a
    case DAE.BINARY(b, DAE.MUL(_), fa)
      guard expHasCref(fa, inExp3) and (not expHasCref(b, inExp3)) and (not expHasCref(inExp2, inExp3))
      equation

        eWithX = Expression.expandFactors(inExp1);
        (factorWithX, factorWithoutX) = List.split1OnTrue(eWithX, expHasCref, inExp3);
        pWithX = makeProductLstSort(factorWithX);
        pWithoutX = makeProductLstSort(factorWithoutX);

        e = Expression.makeDiv(inExp2, pWithoutX);

       then(pWithX, e, true);

    // b*a = rhs  => a = rhs/b solve for a
    case DAE.BINARY(b, DAE.MUL(_), fa)
      guard expHasCref(fa, inExp3) and (not expHasCref(b, inExp3)) and (not expHasCref(inExp2, inExp3))
      equation
        e = Expression.makeDiv(inExp2, b);
       then(fa, e, true);

    // a*b = rhs  => a = rhs/b solve for a
    case DAE.BINARY(fa, DAE.MUL(_), b)
      guard expHasCref(fa, inExp3) and (not expHasCref(b, inExp3)) and (not expHasCref(inExp2, inExp3))
      equation
        e = Expression.makeDiv(inExp2, b);
       then(fa, e, true);

    // f(a)/b = rhs  => f(a) = rhs*b solve for a
    case DAE.BINARY(fa, DAE.DIV(_), b)
      guard expHasCref(fa, inExp3) and (not expHasCref(b, inExp3)) and (not expHasCref(inExp2, inExp3))
      equation
        e = Expression.expMul(inExp2, b);
       then (fa, e, true);

    // g(a)/f(a) = rhs  => rhs*f(a) - g(a) = 0  solve for a
    case DAE.BINARY(ga, DAE.DIV(tp), fa)
      guard expHasCref(fa, inExp3) and expHasCref(ga, inExp3) and (not expHasCref(inExp2, inExp3))
      equation

        e = Expression.expMul(inExp2, fa);
        lhs = Expression.expSub(e, ga);
        e = Expression.makeConstZero(tp);

       then(lhs, e, true);

    else (inExp1, inExp2, false);

  end match;

end preprocessingSolve2;

protected function preprocessingSolve3
"
 helprer function for preprocessingSolve

 (r1)^f(a) = r2 => f(a)  = ln(r2)/ln(r1)
 f(a)^b = 0 => f(a) = 0
 f(a)^n = c => f(a) = c^(1/n)
 abs(x) = 0
 author: Vitalij Ruge
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";

  output DAE.Exp olhs;
  output DAE.Exp orhs;
  output Boolean con "continue";

algorithm
  (olhs, orhs, con) := match(inExp1, inExp2)
      local
       Real r, r1, r2;
       DAE.Exp e1, e2, res;

      // (r1)^f(a) = r2 => f(a)  = ln(r2)/ln(r1)
      case (DAE.BINARY(e1 as DAE.RCONST(r1),DAE.POW(_),e2), DAE.RCONST(r2))
        guard r2 > 0.0 and r1 > 0.0 and (not Expression.isConstOne(e1)) and expHasCref(e2, inExp3)
       equation
         r = log(r2) / log(r1);
         res = DAE.RCONST(r);
       then
         (e2, res, true);

      // f(a)^b = 0 => f(a) = 0
      case (DAE.BINARY(e1,DAE.POW(_),e2), DAE.RCONST(real = 0.0))
        guard expHasCref(e1, inExp3) and (not expHasCref(e2, inExp3))
       then
         (e1, inExp2, true);

      // f(a)^n = c => f(a) = c^(1/n)
      // where n is odd
      case (DAE.BINARY(e1,DAE.POW(_),e2 as DAE.RCONST(r)), _)
        guard (not expHasCref(inExp2, inExp3)) and expHasCref(e1, inExp3) and (1.0 == realMod(r,2.0))
        equation
          res = Expression.makeDiv(DAE.RCONST(1.0),e2);
          res = Expression.expPow(inExp2,res);
       then
         (e1, res, true);

      // sqrt(f(a)) = f(a)^n = c => f(a) = c^(1/n)
      case (DAE.BINARY(e1,DAE.POW(_),DAE.RCONST(0.5)), _)
        guard not expHasCref(inExp2, inExp3) and expHasCref(e1, inExp3)
        equation
          res = Expression.expPow(inExp2,DAE.RCONST(2.0));
       then
         (e1, res, true);

      // abs(x) = 0
      case (DAE.CALL(path = Absyn.IDENT(name = "abs"),expLst = {e1}), DAE.RCONST(0.0))
        then (e1,inExp2,true);

      // sign(x) = 0
      case (DAE.CALL(path = Absyn.IDENT(name = "sign"),expLst = {e1}), DAE.RCONST(0.0))
        then (e1,inExp2,true);


      else (inExp1, inExp2, false);

  end match;


end preprocessingSolve3;


protected function preprocessingSolve4

"
 helprer function for preprocessingSolve

 e.g.
  sqrt(f(x)) - sqrt(g(x))) = 0 = f(x) - g(x)
  exp(f(x)) - exp(g(x))) = 0 = f(x) - g(x)

 author: Vitalij Ruge
"

  input DAE.Exp inExp1;
  input DAE.Exp inExp2;
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  output DAE.Exp oExp1;
  output DAE.Exp oExp2;
  output Boolean newX;

algorithm

  (oExp1, oExp2, newX) := match(inExp1, inExp2, inExp3)
          local
          String s1,s2;
          DAE.Operator op;
          DAE.Exp e1,e2,e3,e4, e, e_1, e_2;
          DAE.Type tp;

          // exp(f(x)) - exp(g(x)) = 0
          case(DAE.BINARY(DAE.CALL(path = Absyn.IDENT("exp"), expLst={e1}), DAE.SUB(_),
                          DAE.CALL(path = Absyn.IDENT("exp"), expLst={e2})),DAE.RCONST(0.0),_)
          then (e1, e2, true);

          // log(f(x)) - log(g(x)) = 0
          case(DAE.BINARY(DAE.CALL(path = Absyn.IDENT("log"), expLst={e1}), DAE.SUB(_),
                          DAE.CALL(path = Absyn.IDENT("log"), expLst={e2})),DAE.RCONST(0.0),_)
          then (e1, e2, true);

          // log10(f(x)) - log10(g(x)) = 0
          case(DAE.BINARY(DAE.CALL(path = Absyn.IDENT("log10"), expLst={e1}), DAE.SUB(_),
                          DAE.CALL(path = Absyn.IDENT("log10"), expLst={e2})),DAE.RCONST(0.0),_)
          then (e1, e2, true);

          // sinh(f(x)) - sinh(g(x)) = 0
          case(DAE.BINARY(DAE.CALL(path = Absyn.IDENT("sinh"), expLst={e1}), DAE.SUB(_),
                          DAE.CALL(path = Absyn.IDENT("sinh"), expLst={e2})),DAE.RCONST(0.0),_)
          then (e1, e2, true);

          // tanh(f(x)) - tanh(g(x)) = 0
          case(DAE.BINARY(DAE.CALL(path = Absyn.IDENT("tanh"), expLst={e1}), DAE.SUB(_),
                          DAE.CALL(path = Absyn.IDENT("tanh"), expLst={e2})),DAE.RCONST(0.0),_)
          then (e1, e2, true);

          // sqrt(f(x)) - sqrt(g(x)) = 0
          case(DAE.BINARY(DAE.CALL(path = Absyn.IDENT("sqrt"), expLst={e1}), DAE.SUB(_),
                          DAE.CALL(path = Absyn.IDENT("sqrt"), expLst={e2})),DAE.RCONST(0.0),_)
          then (e1, e2, true);

          // sinh(f(x)) - cosh(g(x)) = 0
          case(DAE.BINARY(DAE.CALL(path = Absyn.IDENT("sinh"), expLst={e1}), DAE.SUB(_),
                          DAE.CALL(path = Absyn.IDENT("cosh"), expLst={e2})),DAE.RCONST(0.0),_)
          guard Expression.expEqual(e1,e2)
          then (e1, inExp2, true);

          case(DAE.BINARY(DAE.CALL(path = Absyn.IDENT("cosh"), expLst={e1}), DAE.SUB(_),
                          DAE.CALL(path = Absyn.IDENT("sinh"), expLst={e2})),DAE.RCONST(0.0),_)
          guard Expression.expEqual(e1,e2)
          then (e1, inExp2, true);


         // y*sinh(x) - z*cosh(x) = 0
          case(DAE.BINARY(DAE.BINARY(e3,DAE.MUL(),DAE.CALL(path = Absyn.IDENT("sinh"), expLst={e1})), DAE.SUB(tp),
                          DAE.BINARY(e4,DAE.MUL(),DAE.CALL(path = Absyn.IDENT("cosh"), expLst={e2}))),DAE.RCONST(0.0),_)
          guard Expression.expEqual(e1,e2)
          equation
          e = Expression.makePureBuiltinCall("tanh",{e1},tp);
          then (Expression.expMul(e3,e), e4, true);

          case(DAE.BINARY(DAE.BINARY(e4,DAE.MUL(),DAE.CALL(path = Absyn.IDENT("cosh"), expLst={e2})), DAE.SUB(tp),
                          DAE.BINARY(e3,DAE.MUL(),DAE.CALL(path = Absyn.IDENT("sinh"), expLst={e1}))),DAE.RCONST(0.0),_)
          guard Expression.expEqual(e1,e2)
          equation
          e = Expression.makePureBuiltinCall("tanh",{e1},tp);
          then (Expression.expMul(e3,e), e4, true);



          // sqrt(x) - x = 0 -> x = x^2
          case(DAE.BINARY(DAE.CALL(path = Absyn.IDENT("sqrt"), expLst={e1}), DAE.SUB(_),e2), DAE.RCONST(0.0),_)
          then (e1, Expression.expPow(e2, DAE.RCONST(2.0)), true);

          case(DAE.BINARY(e2, DAE.SUB(_),DAE.CALL(path = Absyn.IDENT("sqrt"), expLst={e1})), DAE.RCONST(0.0),_)
          then (e1, Expression.expPow(e2, DAE.RCONST(2.0)), true);

          // f(x)^n - g(x)^n = 0 -> (f(x)/g(x))^n = 1
          case(DAE.BINARY(DAE.BINARY(e1, DAE.POW(), e2), DAE.SUB(tp), DAE.BINARY(e3, DAE.POW(), e4)), DAE.RCONST(0.0),_)
            guard Expression.expEqual(e2,e4) and expHasCref(e1,inExp3) and expHasCref(e3,inExp3)
          equation
            e = Expression.expPow(Expression.makeDiv(e1,e3),e2);
            (e_1, e_2, _) = preprocessingSolve3(e, Expression.makeConstOne(tp), inExp3);
          then (e_1, e_2, true);

          else (inExp1, inExp2, false);

    end match;


end preprocessingSolve4;

protected function expAddX
"
 helprer function for preprocessingSolve

 if(y,g(x),h(x)) + x => if(y, g(x) + x, h(x) + x)
 a*f(x) + b*f(x) = (a+b)*f(x)
 author: Vitalij Ruge
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";

  output DAE.Exp ores;

algorithm
 ores := matchcontinue(inExp1, inExp2, inExp3)
   local
     DAE.Exp e, e1, e2, e3, e4, res;

    case(DAE.IFEXP(e,e1,e2), _,_)
      guard expHasCref(e1, inExp3) and expHasCref(e2, inExp3) and  (not expHasCref(e, inExp3))
     equation
         e3 = expAddX(inExp2, e1, inExp3);
         e4 = expAddX(inExp2, e2, inExp3);

         res = DAE.IFEXP(e, e3, e4);
     then res;

    case(_, DAE.IFEXP(e,e1,e2), _)
      guard expHasCref(e1, inExp3) and expHasCref(e2, inExp3) and (not expHasCref(e, inExp3))
     equation
         e3 = expAddX(inExp1, e1, inExp3);
         e4 = expAddX(inExp1, e2, inExp3);

         res = DAE.IFEXP(e, e3, e4);
     then res;

     else
      equation
       res = expAddX2(inExp1, inExp2, inExp3);
      then res;

 end matchcontinue;

end expAddX;

protected function expAddX2
"
 helprer function for preprocessingSolve
 a*f(x) + b*f(x) = (a+b)*f(x)
 author: Vitalij Ruge
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";

  output DAE.Exp ores;

protected
  list<DAE.Exp> f1, f2;
  DAE.Exp e0,e1,e2;
  DAE.Boolean neg;
  list<DAE.Exp> factorWithX1, factorWithoutX1,  factorWithX2, factorWithoutX2;
  DAE.Exp pWithX1, pWithoutX1, pWithX2, pWithoutX2;

algorithm
  (e0, e1, neg) := match(inExp1)
                   local DAE.Exp ee1, ee2;
                   case(DAE.BINARY(ee1,DAE.ADD(),ee2))
                    then(ee1, ee2, false);
                   case(DAE.BINARY(ee1,DAE.SUB(),ee2))
                    then(ee1, ee2, true);
                   else
                    then(DAE.RCONST(0.0), inExp1, false);
                   end match;

  f1 := Expression.expandFactors(e1);
  (factorWithX1, factorWithoutX1) := List.split1OnTrue(f1, expHasCref, inExp3);
  pWithX1 := makeProductLstSort(factorWithX1);
  pWithoutX1 := makeProductLstSort(factorWithoutX1);
  f2 := Expression.expandFactors(inExp2);
  (factorWithX2, factorWithoutX2) := List.split1OnTrue(f2, expHasCref, inExp3);
  (pWithX2,_) := ExpressionSimplify.simplify1(makeProductLstSort(factorWithX2));
  pWithoutX2 := makeProductLstSort(factorWithoutX2);
  //print("\nf1 =");print(ExpressionDump.printExpListStr(f1));
  //print("\nf2 =");print(ExpressionDump.printExpListStr(f2));

  if Expression.expEqual(pWithX2,pWithX1) then
    // e0 + a*x + b*x -> e0 + (a+b)*x
    if not neg then
      ores := Expression.expAdd(pWithoutX1, pWithoutX2);
    else
    // e0 - a*x + b*x -> e0 + (b-a)*x
      ores := Expression.expSub(pWithoutX2, pWithoutX1);
    end if;
    ores := Expression.expMul(ores, pWithX2);
  elseif Expression.expEqual(pWithX2, Expression.negate(pWithX1)) then
    // e0 + a*(-x) + b*x -> e0 + (b-a)*x
    if not neg then
      ores := Expression.expSub(pWithoutX2, pWithoutX1);
    else
    // e0 - a*(-x) + b*x -> e0 + (b-a)*x
      ores := Expression.expAdd(pWithoutX1, pWithoutX2);
    end if;
    ores := Expression.expMul(ores, pWithX2);
  else
    e1 := Expression.expMul(pWithoutX1, pWithX1);
    e2 := Expression.expMul(pWithoutX2, pWithX2);
    ores := Expression.expAdd(e1,e2);
  end if;

  ores := Expression.expAdd(e0,ores);

end expAddX2;

public function collectX
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp3 "DAE.CREF";
  input DAE.Boolean expand = true;
  output DAE.Exp outLhs;
  output DAE.Exp outRhs;
algorithm
 (outLhs, outRhs) := preprocessingSolve5(inExp1, inExp3, expand);
end collectX;

protected function preprocessingSolve5
"
 helprer function for preprocessingSolve
 split and sort with respect to x
 where x = cref

 f(x,y) = {h(y)*g(x,y), k(y)}

 author: Vitalij Ruge
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  input DAE.Boolean expand;
  output DAE.Exp outLhs;
  output DAE.Exp outRhs;

protected
  list<DAE.Exp> lhs, rhs;
  DAE.Exp tmpLhs, e1;
  Boolean b;
  DAE.ComponentRef cr;
algorithm

   if expHasCref(inExp1, inExp3) then

    if expand then
      (cr, b) := Expression.expOrDerCref(inExp3);
      if b then
        (lhs, rhs) := Expression.allTermsForCref(inExp1, cr, Expression.expHasDerCref);
      else
        (lhs, rhs) := Expression.allTermsForCref(inExp1, cr, Expression.expHasCrefNoPreOrStart);
      end if;
    else
      (lhs, rhs) := List.split1OnTrue(Expression.terms(inExp1), expHasCref, inExp3);
    end if;

    // sort
    // a*f(x)*b -> c*f(x)
    outLhs := DAE.RCONST(0.0);
    tmpLhs := DAE.RCONST(0.0);
    for e in lhs loop
      if Expression.isNegativeUnary(e) then
        DAE.UNARY(exp = e1) := e;
        tmpLhs := expAddX(e1, tmpLhs, inExp3); // special add
      else
        outLhs := expAddX(e, outLhs, inExp3); // special add
      end if;
    end for;
    outLhs := expAddX(outLhs, Expression.negate(tmpLhs), inExp3);
    //rhs
    outRhs := Expression.makeSum1(rhs);
    (outRhs,_) := ExpressionSimplify.simplify1(outRhs);
    (outLhs,_) := ExpressionSimplify.simplify1(outLhs);

   else
    outLhs := DAE.RCONST(0.0);
    outRhs := inExp1;
   end if;

end preprocessingSolve5;

protected function unifyFunCalls
"
e.g.
 semiLinear() -> if
 author: Vitalij Ruge
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  output DAE.Exp oExp;
  output Boolean newX;
algorithm
 (oExp,_) := Expression.traverseExpTopDown(inExp1, unifyFunCallsWork, (inExp3));
 newX := Expression.expEqual(oExp, inExp1);
end unifyFunCalls;

protected function unifyFunCallsWork
  input DAE.Exp inExp;
  input DAE.Exp iT;
  output DAE.Exp outExp;
  output Boolean cont;
  output DAE.Exp oT;
 algorithm
   (outExp,cont,oT) := match(inExp, iT)
   local
     DAE.Exp e, e1,e2,e3, X;
     DAE.Type tp;
/*
   case(DAE.CALL(path = Absyn.IDENT(name = "smooth"), expLst = {_, e}),X)
     guard expHasCref(e, X)
     then (e, true, iT);

   case(DAE.CALL(path = Absyn.IDENT(name = "noEvent"), expLst = {e}),X)
     guard expHasCref(e, X)
     then (e, true, iT);
*/
   case(DAE.CALL(path = Absyn.IDENT(name = "semiLinear"),expLst = {e1, e2, e3}),_)
     guard not Expression.isZero(e1)
       equation
       tp = Expression.typeof(e1);
       e = DAE.IFEXP(DAE.RELATION(e1,DAE.GREATEREQ(tp), Expression.makeConstZero(tp),-1,NONE()),Expression.expMul(e1,e2), Expression.expMul(e1,e3));
     then (e,true, iT);

   // df_der(x) = (x-old(x))/dt
   case(DAE.CALL(path = Absyn.IDENT(name = "$_DF$DER"),expLst = {e1}),X)
     guard expHasCref(e1, X)
     equation
      tp = Expression.typeof(e1);
      e2 = Expression.crefExp(ComponentReference.makeCrefIdent(BackendDAE.symSolverDT, DAE.T_REAL_DEFAULT, {}));
      e3 = Expression.makePureBuiltinCall("pre", {e1}, tp);
      e3 = Expression.expSub(e1,e3);
      e = Expression.expDiv(e3,e2);
     then (e,true, iT);

   else (inExp, true, iT);
   end match;

end unifyFunCallsWork;


protected function solveFunCalls
"
  - inline modelica functions
  - TODO: support annotation inverse
 author: Vitalij Ruge
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  input Option<DAE.FunctionTree> functions;
  output DAE.Exp x;
  output Boolean con;
algorithm
 (x,con) := matchcontinue(functions, inExp1)
                  local DAE.Exp funX; Boolean b;
                  case(_,_)
                  equation
                    (funX,_) = Expression.traverseExpTopDown(inExp1, inlineCallX, (inExp3, functions));
                    b = not Expression.expEqual(funX, inExp1);
                  then (funX, b);
                  else (inExp1, false);
                  end matchcontinue;
end solveFunCalls;

protected function removeSimpleCalls
"
 helprer function for preprocessingSolve

 solve e.g.
   exp(x) = y
   log(x) = y
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";

  output DAE.Exp outLhs;
  output DAE.Exp outRhs;
  output Boolean con "continue";
algorithm
  (outLhs, outRhs, con) := match(inExp1, inExp2, inExp3)
                            case(DAE.CALL(),_,_) then removeSimpleCalls2(inExp1, inExp2, inExp3);
                            else (inExp1, inExp2, false);
                           end match;
end removeSimpleCalls;


protected function removeSimpleCalls2
"
 helprer function for preprocessingSolve

 solve e.g.
   exp(x) = y
   log(x) = y
"
  input DAE.Exp inExp1 "lhs";
  input DAE.Exp inExp2 "rhs";
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";

  output DAE.Exp outLhs;
  output DAE.Exp outRhs;
  output Boolean con "continue";
algorithm
  (outLhs, outRhs, con) := matchcontinue (inExp1,inExp2,inExp3)
    local
      DAE.Exp e1, e2, e3;


    //tanh(x) =y -> x = 1/2 * ln((1+y)/(1-y))
    case (DAE.CALL(path = Absyn.IDENT(name = "tanh"),expLst = {e1}),_,_)
       equation
         true = expHasCref(e1, inExp3);
         false = expHasCref(inExp2, inExp3);
         true = not(Expression.isCref(inExp2) or Expression.isConst(inExp2));
         e2 = Expression.expAdd(DAE.RCONST(1.0), inExp2);
         e3 = Expression.expSub(DAE.RCONST(1.0), inExp2);
         e2 = Expression.makeDiv(e2, e3);
         e2 = Expression.makePureBuiltinCall("log",{e2},DAE.T_REAL_DEFAULT);
         e2 = Expression.expMul(DAE.RCONST(0.5), e2);
       then (e1, e2, true);
    // sinh(x) -> ln(y+(sqrt(1+y^2))
    case (DAE.CALL(path = Absyn.IDENT(name = "sinh"),expLst = {e1}),_,_)
      equation
         true = expHasCref(e1, inExp3);
         false = expHasCref(inExp2, inExp3);
         true = not(Expression.isCref(inExp2) or Expression.isConst(inExp2));
         e2 = Expression.expPow(inExp2, DAE.RCONST(2.0));
         e3 = Expression.expAdd(e2,DAE.RCONST(1.0));
         e2 = Expression.makePureBuiltinCall("sqrt",{e3},DAE.T_REAL_DEFAULT);
         e3 = Expression.expAdd(inExp2, e2);
         e2 = Expression.makePureBuiltinCall("log",{e3},DAE.T_REAL_DEFAULT);
      then (e1,e2,true);

    // log10(f(a)) = g(b) => f(a) = 10^(g(b))
    case (DAE.CALL(path = Absyn.IDENT(name = "log10"),expLst = {e1}),_,_)
       equation
         true = expHasCref(e1, inExp3);
         false = expHasCref(inExp2, inExp3);
         e2 = Expression.expPow(DAE.RCONST(10.0), inExp2);
       then (e1, e2, true);
    // log(f(a)) = g(b) => f(a) = exp(g(b))
    case (DAE.CALL(path = Absyn.IDENT(name = "log"),expLst = {e1}),_,_)
       equation
         true = expHasCref(e1, inExp3);
         false = expHasCref(inExp2, inExp3);
         e2 = Expression.makePureBuiltinCall("exp",{inExp2},DAE.T_REAL_DEFAULT);
       then (e1, e2, true);
    // exp(f(a)) = g(b) => f(a) = log(g(b))
    case (DAE.CALL(path = Absyn.IDENT(name = "exp"),expLst = {e1}),_,_)
       equation
         true = expHasCref(e1, inExp3);
         false = expHasCref(inExp2, inExp3);
         e2 = Expression.makePureBuiltinCall("log",{inExp2},DAE.T_REAL_DEFAULT);
       then (e1, e2, true);
    // sqrt(f(a)) = g(b) => f(a) = (g(b))^2
    case (DAE.CALL(path = Absyn.IDENT(name = "sqrt"),expLst = {e1}),_,_)
       equation
         true = expHasCref(e1, inExp3);
         false = expHasCref(inExp2, inExp3);
         e2 = DAE.RCONST(2.0);
         e2 = Expression.expPow(inExp2,e2);
       then (e1, e2, true);
    // semiLinear(0, a, b) = 0 => a = b // rule 1
    case (DAE.CALL(path = Absyn.IDENT(name = "semiLinear"),expLst = {DAE.RCONST(real = 0.0), e1, e2}),DAE.RCONST(real = 0.0),_)
       then (e1,e2,true);
    // smooth(i,f(a)) = rhs -> f(a) = rhs
    //case (DAE.CALL(path = Absyn.IDENT(name = "smooth"),expLst = {_, e2}),_,_)
    //   then (e2, inExp2, true);
    // noEvent(f(a)) = rhs -> f(a) = rhs
    //case (DAE.CALL(path = Absyn.IDENT(name = "noEvent"),expLst = {e2}),_,_)
    //   then (e2, inExp2, true);

    else (inExp1, inExp2, false);
  end matchcontinue;
end removeSimpleCalls2;

protected function inlineCallX
"
inline function call if depends on X where X is cref or der(cref)
DAE.Exp inExp2 DAE.CREF or 'der(DAE.CREF())'
author: vitalij
"
  input DAE.Exp inExp;
  input tuple<DAE.Exp, Option<DAE.FunctionTree>> iT;
  output DAE.Exp outExp;
  output Boolean cont;
  output tuple<DAE.Exp, Option<DAE.FunctionTree>> oT;
 algorithm
   (outExp,cont,oT) := matchcontinue(inExp, iT)
   local
     DAE.Exp e, X;
     DAE.ComponentRef cr;
     Option<DAE.FunctionTree> functions;
     Boolean b;

   case(DAE.CALL(),(X, functions))
     guard expHasCref(inExp, X)
     equation
       //print("\nfIn: ");print(ExpressionDump.printExpStr(inExp));
       (e,_,b) = Inline.forceInlineExp(inExp,(functions,{DAE.NORM_INLINE(),DAE.DEFAULT_INLINE()}),DAE.emptyElementSource);
       //print("\nfOut: ");print(ExpressionDump.printExpStr(e));
     then (e, not b, iT);
   else (inExp, true, iT);
   end matchcontinue;
end inlineCallX;

protected function preprocessingSolveTmpVars
"
helper function for solveWork
creat tmp vars if needed!
e.g. for solve abs
"
  input DAE.Exp inExp1;
  input DAE.Exp inExp2;
  input DAE.Exp inExp3;
  input Option<DAE.Exp> optCond "condition from an if expression";
  input Integer uniqueEqIndex "offset for tmp vars";
  input list<BackendDAE.Equation> ieqnForNewVars;
  input list<DAE.ComponentRef> inewVarsCrefs;
  input Integer idepth "depth of tmp var";
  output DAE.Exp x;
  output DAE.Exp y;
  output Boolean new_x;
  output list<BackendDAE.Equation> eqnForNewVars "eqn for tmp vars";
  output list<DAE.ComponentRef> newVarsCrefs;
  output Integer odepth;
algorithm
  (x, y, new_x, eqnForNewVars, newVarsCrefs, odepth) := matchcontinue inExp1
    local
      DAE.Exp arg, e1, e_1, e, e2, exP, lhs, e3, e4, e5, e6, rhs, a1,x1, a2,x2, ee1, ee2;
      tuple<DAE.Exp, DAE.Exp> a, c;
      list<DAE.Exp> z1, z2, z3, z4;
      DAE.ComponentRef cr;
      DAE.Type tp;
      BackendDAE.Equation eqn;
      list<BackendDAE.Equation> eqnForNewVars_;
      list<DAE.ComponentRef> newVarsCrefs_;
      DAE.Operator op1, op2;
      String name;

    // try to invert a function call
    // f(x) = y -> x = f^(-1)(y)
    case DAE.CALL(path = Absyn.IDENT(name = name),expLst = {arg})
      guard expHasCref(arg, inExp3) and not expHasCref(inExp2, inExp3)
      algorithm
        (y, new_x, eqnForNewVars_, newVarsCrefs_, odepth) := preprocessingSolveFunctionCall(name, arg, inExp2, inExp3, optCond, uniqueEqIndex, idepth);

        if listEmpty(eqnForNewVars_) then
          eqnForNewVars_ := ieqnForNewVars;
        else
          eqnForNewVars_ := match optCond
            local
              DAE.Exp cond;
            case SOME(cond)
            then BackendDAE.IF_EQUATION({cond}, {eqnForNewVars_}, {}, DAE.emptyElementSource, BackendDAE.EQ_ATTR_DEFAULT_UNKNOWN) :: ieqnForNewVars;
            else listAppend(eqnForNewVars_, ieqnForNewVars);
          end match;
        end if;
      then (if new_x then arg else inExp1, y, new_x, eqnForNewVars_, listAppend(newVarsCrefs_, inewVarsCrefs), odepth);

    // x^n = y -> x = y^(1/n)
    case DAE.BINARY(e1, DAE.POW(tp), e2)
      guard expHasCref(e1, inExp3) and not expHasCref(e2, inExp3)
      algorithm
        tp := Expression.typeof(e1);
        exP := makeInitialGuess(tp, inExp3, e1);
        // exP := makeInitialGuess(tp, inExp3, inExp2);
        (exP, eqnForNewVars_, newVarsCrefs_) := makeTmpEqnAndCrefFromExp(exP, tp, "X$ABS", uniqueEqIndex, idepth, ieqnForNewVars, inewVarsCrefs, false);
        e_1 := Expression.makePureBuiltinCall("$_signNoNull", {exP}, tp); // sign

        lhs := Expression.expPow(inExp2, Expression.inverseFactors(e2)); // y^(1/n)
        lhs := Expression.makePureBuiltinCall("abs", {lhs}, tp); // abs(y^(1/n))
        //lhs := Expression.makePureBuiltinCall("abs", {inExp2}, tp); // abs(y)
        //lhs := Expression.expPow(lhs, Expression.inverseFactors(e2)); // abs(y)^(1/n)
        lhs := Expression.expMul(e_1, lhs); // sign*abs(y^(1/n))
      then(e1, lhs, true, eqnForNewVars_, newVarsCrefs_, idepth + 1);

    //QE
    // a*x^n + b*x^m = c
    // a*x^n - b*x^m = c
    case DAE.BINARY(ee1, op1, ee2)
      guard(Expression.isAddOrSub(op1))
      algorithm
        (z1, z2) := List.split1OnTrue(Expression.factors(ee1), expHasCref, inExp3);
        (z3, z4) := List.split1OnTrue(Expression.factors(ee2), expHasCref, inExp3);

        x1 := makeProductLstSort(z1);
        a1 := makeProductLstSort(z2);

        x2 := makeProductLstSort(z3);
        a2 := if Expression.isAdd(op1) then makeProductLstSort(z4) else Expression.negate(makeProductLstSort(z4));
        a := simplifyBinaryMulCoeff(x1);
        c := simplifyBinaryMulCoeff(x2);
        (e2, e3) := a;
        (e5, e6) := c;
        (lhs, rhs, eqnForNewVars_, newVarsCrefs_) := solveQE(a1,e2,e3,a2,e5,e6,inExp2,inExp3,ieqnForNewVars,inewVarsCrefs,uniqueEqIndex,idepth);
      then(lhs, rhs, true, eqnForNewVars_, newVarsCrefs_, idepth + 1);

    else (inExp1, inExp2, false, ieqnForNewVars, inewVarsCrefs, idepth);
  end matchcontinue;
end preprocessingSolveTmpVars;

protected function preprocessingSolveFunctionCall
  input String name "of the function";
  input DAE.Exp arg "ument of the function";
  input DAE.Exp rhs;
  input DAE.Exp inExp3 "solve for this";
  input Option<DAE.Exp> optCond "condition from an if expression";
  input Integer uniqueEqIndex "offset for tmp vars";
  input Integer idepth "depth of tmp var";
  output DAE.Exp result;
  output Boolean new_x;
  output list<BackendDAE.Equation> newEqns "eqns for tmp vars";
  output list<DAE.ComponentRef> newVars "tmp vars";
  output Integer odepth;
algorithm
  (result, new_x, newEqns, newVars, odepth) := match name
    local
      DAE.Exp y, exP, sgn, inv, pi, e1, e2, k1, k2, x1, x2;
      DAE.Type tp;
      list<BackendDAE.Equation> eqns;
      list<DAE.ComponentRef> vars;
      BackendDAE.Equation ass "assertion for finding domain violations";

    // tanh(x) -> 0.5 * ln((1 + y)/(1 - y))
    // exists for y in (-1, 1)
    // unique
    case "tanh" algorithm
      tp := Expression.typeof(rhs);
      (y, eqns, vars) := makeTmpEqnAndCrefFromExp(rhs, tp, "Y$TANH", uniqueEqIndex, idepth, {}, {}, false);
      e1 := Expression.expAdd(DAE.RCONST(1.0), y); // 1 + y
      e2 := Expression.expSub(DAE.RCONST(1.0), y); // 1 - y
      e1 := Expression.makeDiv(e1, e2); // (1 + y)/(1 - y)
      e1 := Expression.makePureBuiltinCall("log", {e1}, tp); // ln((1 + y)/(1 - y))
      inv := Expression.expMul(DAE.RCONST(0.5), e1); // 0.5 * ln((1 + y)/(1 - y))
      ass := makeDomainAssert(name, rhs, SOME((-1.0, false)), SOME((1.0, false))); // y in (-1, 1)
    then (inv, true, ass :: eqns, vars, idepth + 1);

    // sinh(x) -> ln(y + sqrt(y^2 + 1))
    // exixts always
    // unique
    case "sinh" algorithm
      tp := Expression.typeof(rhs);
      (y, eqns, vars) := makeTmpEqnAndCrefFromExp(rhs, tp, "Y$SINH", uniqueEqIndex, idepth, {}, {}, false);
      e1 := Expression.expPow(y, DAE.RCONST(2.0)); // y^2
      e1 := Expression.expAdd(e1, DAE.RCONST(1.0)); // y^2 + 1
      e1 := Expression.makePureBuiltinCall("sqrt", {e1}, tp); // sqrt(y^2 + 1)
      e1 := Expression.expAdd(y, e1); // y + sqrt(y^2 + 1)
      e1 := Expression.makePureBuiltinCall("log", {e1}, tp); // ln(y + sqrt(y^2 + 1))
    then (e1, true, eqns, vars, idepth + 1);

    // cosh(x) -> ln(y + sign*sqrt(y^2 - 1))
    // exists for y in [1, inf)
    // two values, sign is -1 or 1
    case "cosh" algorithm
      tp := Expression.typeof(rhs);
      (y, eqns, vars) := makeTmpEqnAndCrefFromExp(rhs, tp, "Y$COSH", uniqueEqIndex, idepth, {}, {}, false);

      exP := makeInitialGuess(tp, inExp3, arg);
      (exP, eqns, vars) := makeTmpEqnAndCrefFromExp(exP, tp, "SIGN$COSH", uniqueEqIndex, idepth, eqns, vars, false);
      sgn := Expression.makePureBuiltinCall("$_signNoNull", {exP}, tp); // sign

      e1 := Expression.expPow(y, DAE.RCONST(2.0)); // y^2
      e1 := Expression.expSub(e1, DAE.RCONST(1.0)); // y^2 - 1
      e1 := Expression.makePureBuiltinCall("sqrt", {e1}, tp); // sqrt(y^2 - 1)
      e1 := Expression.expMul(sgn, e1); // sign*sqrt(y^2 - 1)
      e1 := Expression.expAdd(y, e1); // y + sign*sqrt(y^2 - 1)
      e1 := Expression.makePureBuiltinCall("log", {e1}, tp); // ln(y + sign*sqrt(y^2 - 1))

      ass := makeDomainAssert(name, rhs, SOME((1.0, true)), NONE()); // y in [1, inf)
    then (e1, true, ass :: eqns, vars, idepth + 1);

    // cos(x) -> sign*acos(y) + 2*pi*k
    // exists for y in [-1, 1]
    // infinitely many values, k is integer, sign is -1 or 1
    case "cos" algorithm
      tp := Expression.typeof(rhs);
      (y, eqns, vars) := makeTmpEqnAndCrefFromExp(rhs, tp, "Y$COS", uniqueEqIndex, idepth, {}, {}, false);

      inv := Expression.makePureBuiltinCall("acos", {y}, tp);
      (inv, eqns, vars) := makeTmpEqnAndCrefFromExp(inv, tp, "INV$COS", uniqueEqIndex, idepth, eqns, vars, false);

      exP := makeInitialGuess(tp, inExp3, arg);
      (exP, eqns, vars) := makeTmpEqnAndCrefFromExp(exP, tp, "PREX$COS", uniqueEqIndex, idepth, eqns, vars, false);

      k1 := helpInvCos(inv, exP, tp, true);
      k2 := helpInvCos(inv, exP, tp, false);
      (k1, eqns, vars) := makeTmpEqnAndCrefFromExp(k1, tp, "k1$COS", uniqueEqIndex, idepth, eqns, vars, false);
      (k2, eqns, vars) := makeTmpEqnAndCrefFromExp(k2, tp, "k2$COS", uniqueEqIndex, idepth, eqns, vars, false);

      x1 := helpInvCos2(k1, inv, tp, true);
      x2 := helpInvCos2(k2, inv, tp, false);
      (x1, eqns, vars) := makeTmpEqnAndCrefFromExp(x1, tp, "x1$COS", uniqueEqIndex, idepth, eqns, vars, false);
      (x2, eqns, vars) := makeTmpEqnAndCrefFromExp(x2, tp, "x2$COS", uniqueEqIndex, idepth, eqns, vars, false);
      e1 := helpInvCos3(x1, x2, exP, tp);

      ass := makeDomainAssert(name, rhs, SOME((-1.0, true)), SOME((1.0, true))); // y in [-1, 1]
    then (e1, true, ass :: eqns, vars, idepth + 1);

    // sin(x) -> (-1)^k * asin(y) + k*pi
    // exists for y in [-1, 1]
    // infinitely many values, k is integer
    case "sin" algorithm
      tp := Expression.typeof(rhs);
      (y, eqns, vars) := makeTmpEqnAndCrefFromExp(rhs, tp, "Y$SIN", uniqueEqIndex, idepth, {}, {}, false);

      inv := Expression.makePureBuiltinCall("asin", {y}, tp);
      (inv, eqns, vars) := makeTmpEqnAndCrefFromExp(inv, tp, "INV$SIN", uniqueEqIndex, idepth, eqns, vars, false);

      exP := makeInitialGuess(tp,inExp3,arg);
      (exP, eqns, vars) := makeTmpEqnAndCrefFromExp(exP, tp, "PREX$SIN", uniqueEqIndex, idepth, eqns, vars, false);

      k1 := helpInvSin(inv, arg, tp, true);
      k2 := helpInvSin(inv, arg, tp, false);
      (k1, eqns, vars) := makeTmpEqnAndCrefFromExp(k1, tp, "k1$SIN", uniqueEqIndex, idepth, eqns, vars, false);
      (k2, eqns, vars) := makeTmpEqnAndCrefFromExp(k2, tp, "k2$SIN", uniqueEqIndex, idepth, eqns, vars, false);

      x1 := helpInvSin2(k1, inv, tp, true);
      x2 := helpInvSin2(k2, inv, tp, false);
      (x1, eqns, vars) := makeTmpEqnAndCrefFromExp(x1, tp, "x1$SIN", uniqueEqIndex, idepth, eqns, vars, false);
      (x2, eqns, vars) := makeTmpEqnAndCrefFromExp(x2, tp, "x2$SIN", uniqueEqIndex, idepth, eqns, vars, false);

      e1 := helpInvCos3(x1, x2, exP, tp);

      ass := makeDomainAssert(name, rhs, SOME((-1.0, true)), SOME((1.0, true))); // y in [-1, 1]
    then (e1, true, ass :: eqns, vars, idepth + 1);

    // tan(x) -> atan(y) + k*pi
    // exists always
    // infinitely many values, k is integer
    case "tan" algorithm
      tp := Expression.typeof(rhs);
      (y, eqns, vars) := makeTmpEqnAndCrefFromExp(rhs, tp, "Y$TAN", uniqueEqIndex, idepth, {}, {}, false);

      inv := Expression.makePureBuiltinCall("atan", {y}, tp);
      (inv, eqns, vars) := makeTmpEqnAndCrefFromExp(inv, tp, "INV$TAN", uniqueEqIndex, idepth, eqns, vars, false);

      exP := makeInitialGuess(tp, inExp3, arg);
      (exP, eqns, vars) := makeTmpEqnAndCrefFromExp(exP, tp, "PREX$TAN", uniqueEqIndex, idepth, eqns, vars, false);

      k1 := Expression.expSub(exP, inv); // pre(x) - atan(y)
      k1 := Expression.makeDiv(k1, DAE.PI); // (pre(x) - atan(y))/pi
      k1 := Expression.makePureBuiltinCall("$_round", {k1}, tp); // k = round((pre(x) - atan(y))/pi)
      e1 := Expression.expMul(k1, DAE.PI); // k*pi
      e1 := Expression.expAdd(inv, e1); // atan(y) + pi*k
    then (e1, true, eqns, vars, idepth + 1);

    // abs(x) -> sign*y
    // exists for y in [0, inf)
    // two values, sign is -1 or 1
    case "abs" algorithm
      tp := Expression.typeof(arg);
      exP := makeInitialGuess(tp, inExp3, arg);
      (exP, eqns, vars) := makeTmpEqnAndCrefFromExp(exP, tp, "SIGN$ABS", uniqueEqIndex, idepth, {}, {}, false);
      sgn := Expression.makePureBuiltinCall("$_signNoNull", {exP}, tp); // sign
      e1 := Expression.expMul(sgn, rhs); // sign*y
      ass := makeDomainAssert(name, rhs, SOME((0.0, true)), NONE()); // y in [0, inf)
    then (e1, true, ass ::eqns, vars, idepth + 1);

    // sqrt(x) -> y^2
    // exists for y in [0, inf)
    // unique
    case "sqrt" algorithm
      inv := Expression.expPow(rhs, DAE.RCONST(2.0)); // y^2
      ass := makeDomainAssert(name, rhs, SOME((0.0, true)), NONE()); // y in [0, inf)
    then (inv, true, {ass}, {}, idepth + 1);

    // asin(x) -> sin(y)
    // exists for y in [-pi/2, pi/2]
    // unique
    case "asin" algorithm
      tp := Expression.typeof(rhs);
      inv := Expression.makePureBuiltinCall("sin", {rhs}, tp); // sin(y)
      ass := makeDomainAssert(name, rhs, SOME((-0.5*Expression.toReal(DAE.PI), true)), SOME((0.5*Expression.toReal(DAE.PI), true))); // y in [-pi/2, pi/2]
    then (inv, true, {ass}, {}, idepth + 1);

    // acos(x) -> cos(y)
    // exists for y in [0, pi]
    // unique
    case "acos" algorithm
      tp := Expression.typeof(rhs);
      inv := Expression.makePureBuiltinCall("cos", {rhs}, tp); // cos(y)
      ass := makeDomainAssert(name, rhs, SOME((0.0, true)), SOME((Expression.toReal(DAE.PI), true))); // y in [0, pi]
    then (inv, true, {ass}, {}, idepth + 1);

    // atan(x) -> tan(y)
    // exists for y in [-pi/2, pi/2]
    // unique
    case "atan" algorithm
      tp := Expression.typeof(rhs);
      inv := Expression.makePureBuiltinCall("tan", {rhs}, tp); // tan(y)
      ass := makeDomainAssert(name, rhs, SOME((-0.5*Expression.toReal(DAE.PI), true)), SOME((0.5*Expression.toReal(DAE.PI), true))); // y in [-pi/2, pi/2]
    then (inv, true, {ass}, {}, idepth + 1);

    // exp(x) -> log(y)
    // exists for y in (0, inf)
    // unique
    case "exp" algorithm
      tp := Expression.typeof(rhs);
      inv := Expression.makePureBuiltinCall("log", {rhs}, tp); // log(y)
      ass := makeDomainAssert(name, rhs, SOME((0.0, false)), NONE()); // y in (0, inf)
    then (inv, true, {ass}, {}, idepth + 1);

    // log(x) -> exp(y)
    // exists always
    // unique
    case "log" algorithm
      tp := Expression.typeof(rhs);
      inv := Expression.makePureBuiltinCall("exp", {rhs}, tp); // exp(y)
    then (inv, true, {}, {}, idepth + 1);

    // log10(x) -> 10^y
    // exists always
    // unique
    case "log10" algorithm
      inv := Expression.expPow(DAE.RCONST(10.0), rhs); // 10^y
    then (inv, true, {}, {}, idepth + 1);

    // sign(x) is not invertible
    case "sign" then (rhs, false, {}, {}, idepth);

    // $_DF$DER(x) = y  ->  (x - pre(x))/dt = y  ->  x = y*dt + pre(x)
    case "$_DF$DER" algorithm
      e1 := Expression.crefExp(ComponentReference.makeCrefIdent(BackendDAE.symSolverDT, DAE.T_REAL_DEFAULT, {})); // dt
      exP := Expression.makePureBuiltinCall("pre", {arg}, Expression.typeof(arg)); // pre(x)
      e1 := Expression.expAdd(Expression.expMul(rhs, e1), exP); // y*dt + pre(x)
    then(e1, true, {}, {}, idepth + 1);

    // don't know inverse of this function
    else (rhs, false, {}, {}, idepth);
  end match;
end preprocessingSolveFunctionCall;

protected function simplifyBinaryMulCoeff
"generalization of ExpressionSimplify.simplifyBinaryMulCoeff2"
  input DAE.Exp inExp;
  output tuple<DAE.Exp, DAE.Exp> outRes;
algorithm
  outRes := match(inExp)
    local
      DAE.Exp e,e1,e2;
      DAE.Exp coeff;

    case ((e as DAE.CREF()))
      then ((e, DAE.RCONST(1.0)));

    case (DAE.BINARY(exp1 = e1,operator = DAE.POW(),exp2 = DAE.UNARY(operator = DAE.UMINUS(), exp = coeff)))
      then
        ((e1, Expression.negate(coeff)));

    case (DAE.BINARY(exp1 = e1,operator = DAE.POW(),exp2 = coeff))
      then ((e1,coeff));

    case (DAE.BINARY(exp1 = e1,operator = DAE.MUL(),exp2 = e2))
    guard(Expression.expEqual(e1, e2))
      then
        ((e1, DAE.RCONST(2.0)));

    case(DAE.BINARY(e1, DAE.DIV(), e2))
    guard(Expression.isOne(e1))
      then(e2, DAE.RCONST(-1.0));

    case(DAE.CALL(path=Absyn.IDENT("sqrt"),expLst={e}))
      then ((e,DAE.RCONST(0.5)));

    else ((inExp,DAE.RCONST(1.0)));

  end match;
end simplifyBinaryMulCoeff;

protected function solveQE
"
solve Quadratic equation with respect to inExp3
IN: a,x,n,b,y,m
where solve a*x^n + b*y^m = inExp2 with 2*m = n or 2*n = m and y = x

author: Vitalij Ruge
"
 input DAE.Exp e1,e2,e3,e4,e5,e6;
 input DAE.Exp inExp2;
 input DAE.Exp inExp3;

 input list<BackendDAE.Equation> ieqnForNewVars "eqn for tmp vars";
 input list<DAE.ComponentRef> inewVarsCrefs "cref for tmp vars";
 input Integer uniqueEqIndex, idepth "need for tmp vars";

 output DAE.Exp rhs, lhs;
 output list<BackendDAE.Equation> eqnForNewVars;
 output list<DAE.ComponentRef> newVarsCrefs;

protected
  DAE.Exp e, e7, con, invExp, x1, x2, x, exP;
  DAE.Exp a,b,c, n, sgnb, b2, ac, sExp1, sExp2;
  DAE.ComponentRef cr;
  DAE.Type tp;
  BackendDAE.Equation eqn;
  Boolean b1, b3;
algorithm
    false := Expression.isZero(e1) and Expression.isZero(e2);
    true := Expression.expEqual(e2,e5);
    b1 := Expression.expEqual(e3, Expression.expMul(DAE.RCONST(2.0),e6));
    b3 := Expression.expEqual(e6, Expression.expMul(DAE.RCONST(2.0),e3));

    true := b1 or b3;
    false := expHasCref(e1, inExp3);
    true := expHasCref(e2, inExp3);
    false := expHasCref(e3, inExp3);
    false := expHasCref(e4, inExp3);
    true := expHasCref(e5, inExp3);
    false := expHasCref(e6, inExp3);
    false := expHasCref(inExp2, inExp3);


    a := if b1 then e1 else e4;
    b := if b1 then e4 else e1;
    c := Expression.negate(inExp2);
    n := if b1 then e6 else e3;

    tp := Expression.typeof(a);
    (a, eqnForNewVars, newVarsCrefs) := makeTmpEqnAndCrefFromExp(a, tp, "a$QE", uniqueEqIndex, idepth, ieqnForNewVars, inewVarsCrefs, false);
    con := DAE.RELATION(a,DAE.EQUAL(tp),DAE.RCONST(0.0),-1,NONE());

    tp := Expression.typeof(b);
    (b, eqnForNewVars, newVarsCrefs) := makeTmpEqnAndCrefFromExp(b, tp, "b$QE", uniqueEqIndex, idepth, eqnForNewVars, newVarsCrefs, false);
    sgnb := Expression.makePureBuiltinCall("$_signNoNull",{b},tp);
    b2 := Expression.expPow(b, DAE.RCONST(2.0));
    (b2, eqnForNewVars, newVarsCrefs) := makeTmpEqnAndCrefFromExp(b2, tp, "bPow2$QE", uniqueEqIndex, idepth, eqnForNewVars, newVarsCrefs, false);

    tp := Expression.typeof(c);
    (c, eqnForNewVars, newVarsCrefs) := makeTmpEqnAndCrefFromExp(c, tp, "c$QE", uniqueEqIndex, idepth, eqnForNewVars, newVarsCrefs, false);
    ac := Expression.expMul(a,c);
    ac := Expression.expMul(DAE.RCONST(4.0),ac);

    sExp1 := Expression.expSub(b2,ac);
    sExp2 := Expression.makePureBuiltinCall("sqrt",{sExp1},tp);
    sExp2 := Expression.expMul(sgnb, sExp2);


    a := DAE.IFEXP(con, Expression.makeConstOne(tp), a);
    (a, eqnForNewVars, newVarsCrefs) := makeTmpEqnAndCrefFromExp(a, tp, "a1$QE", uniqueEqIndex, idepth, eqnForNewVars, newVarsCrefs, false);

    x1 := Expression.expAdd(b, sExp2);
    x1 := Expression.makeDiv(x1, a);
    x1 := Expression.expMul(DAE.RCONST(-0.5), x1);
    tp := Expression.typeof(x1);
    x1 := DAE.IFEXP(con, Expression.makeConstOne(tp), x1);
    (x1, eqnForNewVars, newVarsCrefs) := makeTmpEqnAndCrefFromExp(x1, tp, "x1$QE", uniqueEqIndex, idepth, eqnForNewVars, newVarsCrefs, false);

    //Vieta
    x2 := Expression.expMul(a,x1);
    x2 := Expression.makeDiv(c,x2);
    x2 := DAE.IFEXP(con, Expression.makeConstOne(tp), x2);
    x2 := DAE.IFEXP(DAE.RELATION(x1,DAE.EQUAL(tp),DAE.RCONST(0.0),-1,NONE()), DAE.RCONST(0.0), x2);
    (x2, eqnForNewVars, newVarsCrefs) := makeTmpEqnAndCrefFromExp(x2, tp, "x2$QE", uniqueEqIndex, idepth, eqnForNewVars, newVarsCrefs, false);

    tp := Expression.typeof(e2);
    exP := makeInitialGuess(tp,inExp3,e2);
    (exP, eqnForNewVars, newVarsCrefs) := makeTmpEqnAndCrefFromExp(exP, tp, "prex$QE", uniqueEqIndex, idepth, eqnForNewVars, newVarsCrefs, false);

    x := helpInvCos3(x1,x2,exP,tp);
    (x, eqnForNewVars, newVarsCrefs) := makeTmpEqnAndCrefFromExp(x, tp, "x$QE", uniqueEqIndex, idepth, eqnForNewVars, newVarsCrefs, false);

    // a = 0
    e7 :=  Expression.makeDiv(inExp2,b);
    invExp := Expression.inverseFactors(n);
    (invExp, _) :=  ExpressionSimplify.simplify1(invExp);
    e7 := Expression.expPow(e7, invExp);

    // if a==0
    rhs := DAE.IFEXP(con, e7 , x);
    // lhs
    lhs := if b1 then Expression.expPow(e2, e6) else Expression.expPow(e2, e3);

end solveQE;

protected function solveIfExp
"
 solve:
  if(f(y), f(x), g(x) ) = h(y) w.r.t. x
"
  input DAE.Exp inExp1;
  input DAE.Exp inExp2;
  input DAE.Exp inExp3;
  input Option<DAE.Exp> inCond;
  input Option<DAE.FunctionTree> functions;
  input Option<Integer> uniqueEqIndex "offset for tmp vars";
  input Integer idepth;
  input Boolean doInline;
  input Boolean isContinuousIntegration;
  output DAE.Exp outExp;
  output list<DAE.Statement> outAsserts;
  output list<BackendDAE.Equation> eqnForNewVars "eqn for tmp vars";
  output list<DAE.ComponentRef> newVarsCrefs;
  output Integer odepth;
algorithm
  (outExp, outAsserts, eqnForNewVars, newVarsCrefs, odepth) := match inExp1
    local
      DAE.Exp eCond, eThen, eElse, res, lhs, rhs, cond1, cond2;
      list<DAE.Statement> asserts, asserts1, asserts2;
      list<BackendDAE.Equation> eqns, eqns1;
      list<DAE.ComponentRef> var, var1;
      Integer depth;

    //  f(a) if(g(b)) then f1(a) else f2(a) =>
    //  a1 = solve f(a),f1(a) for a
    //  a2 = solve f(a),f2(a) for a
    //  => a = if g(b) then a1 else a2
    case DAE.IFEXP(eCond, eThen, eElse)
      guard
        isContinuousIntegration or not expHasCref(eCond, inExp3)
      algorithm

        // nested if expressions need to combine their conditions
        (cond1, cond2) := match inCond
          local DAE.Exp theCond;
          case SOME(theCond)
            then (DAE.LBINARY(theCond, DAE.AND(Expression.typeof(eCond)), eCond), DAE.LBINARY(theCond, DAE.AND(Expression.typeof(eCond)), Expression.negate(eCond)));
          else (eCond, Expression.negate(eCond));
        end match;

        (lhs, asserts1, eqns, var, depth) := solveWork(eThen, inExp2, inExp3, SOME(cond1), functions, uniqueEqIndex, idepth, doInline, isContinuousIntegration);
        (rhs, _, eqns1, var1, depth) := solveWork(eElse, inExp2, inExp3, SOME(cond2), functions, uniqueEqIndex, depth, doInline, isContinuousIntegration);

        res := DAE.IFEXP(eCond, lhs, rhs);
        asserts := listAppend(asserts1, asserts1);
      then
        (res, asserts, listAppend(eqns1, eqns), listAppend(var1, var), depth);
    else fail();
  end match;
end solveIfExp;

protected function solveLinearSystem
"
 solve linear system with newton step

 ToDo:
  fixed is for ./simulation/modelica/equations/deriveToLog.mos
"
  input DAE.Exp inExp1;
  input DAE.Exp inExp2;
  input DAE.Exp inExp3;
  input Option<DAE.FunctionTree> functions;
  input Integer idepth;
  output DAE.Exp outExp;
  output list<DAE.Statement> outAsserts;
  output list<BackendDAE.Equation> eqnForNewVars = {} "eqn for tmp vars";
  output list<DAE.ComponentRef> newVarsCrefs = {};
  output Integer odepth = idepth;


algorithm
   (outExp,outAsserts) := match(inExp1,inExp2,inExp3)
   local
      DAE.Exp dere,e,z;
      DAE.ComponentRef cr;
      DAE.Exp rhs;
      DAE.Type tp;
      Integer i;

    // cr = (e1-e2)/(der(e1-e2,cr))
    case (_,_,DAE.CREF(componentRef = cr))
      equation
        false = hasOnlyFactors(inExp1,inExp2);
        e = Expression.expSub(inExp1,inExp2);
        (e,_) = ExpressionSimplify.simplify1(e);
        //print("\ne: ");print(ExpressionDump.printExpStr(e));
        dere = Differentiate.differentiateExpSolve(e, cr, functions);
        //print("\nder(e): ");print(ExpressionDump.printExpStr(dere));
        (dere,_) = ExpressionSimplify.simplify(dere);
        false = Expression.isZero(dere);
        false = Expression.expHasCrefNoPreOrStart(dere, cr);
        tp = Expression.typeof(inExp3);
        (z,_) = Expression.makeZeroExpression(Expression.arrayDimension(tp));
        ((e,i)) = Expression.replaceExp(e, inExp3, z);
        // replace at least once, otherwise it's wrong
        if i < 1 then
          fail();
        end if;
        (e,_) = ExpressionSimplify.simplify(e);
        rhs = Expression.negate(Expression.makeDiv(e,dere));
      then
        (rhs,{});

      else fail();
   end match;

end solveLinearSystem;

protected function hasOnlyFactors "help function to solve2, returns true if equation e1 == e2, has either e1 == 0 or e2 == 0 and the expression only contains
factors, e.g. a*b*c = 0. In this case we can not solve the equation"
  input DAE.Exp e1;
  input DAE.Exp e2;
  output Boolean res;
algorithm
  res := matchcontinue(e1,e2)

    // try normal
    case(_,_)
      equation
        true = Expression.isZero(e1);
        // More than two factors
        _::_::_ = Expression.factors(e2);
        //.. and more than two crefs
        _::_::_ = Expression.extractCrefsFromExp(e2);
      then
        true;

    // swapped args
    case(_,_)
      equation
        true = Expression.isZero(e2);
        _::_::_ = Expression.factors(e1);
        _::_::_ = Expression.extractCrefsFromExp(e1);
      then
        true;

    else false;

  end matchcontinue;
end hasOnlyFactors;


protected function expHasCref
"
helper function for solve.
case distinction for
DAE.CREF or 'der(DAE.CREF())'
Expression.expHasCrefNoPreOrStart
or
Expression.expHasDerCref
"
  input DAE.Exp inExp1;
  input DAE.Exp inExp3 "DAE.CREF or 'der(DAE.CREF())'";
  output DAE.Boolean res;

algorithm
  res := match(inExp1, inExp3)
         local DAE.ComponentRef cr;

          case(_, DAE.CREF(componentRef = cr)) then Expression.expHasCrefNoPreOrStart(inExp1, cr);
          case(_, DAE.CALL(path = Absyn.IDENT(name = "der"),expLst = {DAE.CREF(componentRef = cr)})) then Expression.expHasDerCref(inExp1, cr);
          else
           equation
            if Flags.isSet(Flags.FAILTRACE) then
              print("\n-ExpressionSolve.solve failed:");
              print(" with respect to: ");print(ExpressionDump.printExpStr(inExp3));
              print(" not support!");
              print("\n");
            end if;
          then fail();
         end match;

end expHasCref;

protected function makeProductLstSort
 "Takes a list of expressions an makes a product
  expression multiplying all elements in the list.

- a*if(b,c,d) -> if(b,a*c,a*d)

"
  input list<DAE.Exp> inExpLst;
  output DAE.Exp outExp;
protected
  DAE.Type tp;
  list<DAE.Exp> expLstDiv, expLst, expLst2;
  DAE.Exp e, e1, e2;
  DAE.Operator op;
algorithm
  if listEmpty(inExpLst) then
    outExp := DAE.RCONST(1.0);
  return;
  end if;

  tp := Expression.typeof(listHead(inExpLst));

  (expLstDiv, expLst) :=  List.splitOnTrue(inExpLst, Expression.isDivBinary);
  outExp := makeProductLstSort2(expLst, tp);
  if not listEmpty(expLstDiv) then
    expLst2 := {};
    expLst := {};

    for elem in expLstDiv loop
      DAE.BINARY(e1,op,e2) := elem;
      expLst := e1::expLst;
      expLst2 := e2::expLst2;
    end for;

    if not listEmpty(expLst2) then
      e := makeProductLstSort(expLst2);
      if not Expression.isOne(e) then
        outExp := Expression.makeDiv(outExp,e);
      end if;
    end if;

    if not listEmpty(expLst) then
      e := makeProductLstSort(expLst);
      outExp := Expression.expMul(outExp,e);
    end if;

  end if;

end makeProductLstSort;


protected function makeProductLstSort2
  input list<DAE.Exp> inExpLst;
  input DAE.Type tp;
  output DAE.Exp outExp = Expression.makeConstOne(tp);
protected
  list<DAE.Exp> rest;
algorithm
  rest := ExpressionSimplify.simplifyList(inExpLst);
  for elem in rest loop
    if not Expression.isOne(elem) then
    outExp := match(elem)
              local DAE.Exp e1,e2,e3;
              case(DAE.IFEXP(e1,e2,e3))
              then DAE.IFEXP(e1, Expression.expMul(outExp,e2), Expression.expMul(outExp,e3));
              else Expression.expMul(outExp, elem);
              end match;
    end if;
  end for;

end makeProductLstSort2;

protected function makeTmpEqnAndCrefFromExp
  input DAE.Exp iExp;
  input DAE.Type tp;
  input String name;
  input Integer index1, index2;
  input list<BackendDAE.Equation> ieqnForNewVars;
  input list<DAE.ComponentRef> inewVarsCrefs;
  input Boolean need;
  output DAE.Exp oExp;
  output list<BackendDAE.Equation> oeqnForNewVars;
  output list<DAE.ComponentRef> onewVarsCrefs;
protected
  DAE.ComponentRef cr;
  BackendDAE.Equation eqn;
algorithm
  (oExp,_) := ExpressionSimplify.simplify1(iExp);
  if need or not (Expression.isCref(oExp) or Expression.isConst(oExp)) then
    cr := ComponentReference.makeCrefIdent("$TMP$VAR$" + intString(index1) + "$" + intString(index2) + name, tp , {});
    eqn := BackendDAE.SOLVED_EQUATION(cr, oExp, DAE.emptyElementSource, BackendDAE.EQ_ATTR_DEFAULT_UNKNOWN);
    oExp := Expression.crefExp(cr);
    oeqnForNewVars := eqn :: ieqnForNewVars;
    onewVarsCrefs := cr :: inewVarsCrefs;
  else
    oeqnForNewVars := ieqnForNewVars;
    onewVarsCrefs := inewVarsCrefs;
  end if;
end makeTmpEqnAndCrefFromExp;

protected function makeDomainAssert
  input String name "of the function";
  input DAE.Exp rhs "solution of the function";
  input Option<tuple<Real, Boolean>> lowerBound "(value, including?)";
  input Option<tuple<Real, Boolean>> upperBound "(value, including?)";
  output BackendDAE.Equation assEq;
protected
  String msg;
  DAE.Exp cond;
  DAE.Algorithm algo;
  DAE.Type tp = Expression.typeof(rhs);
algorithm
  (msg, cond) := match (lowerBound, upperBound)
    local
      Real lower, upper;
      String str;
      DAE.Exp l, u;

    // range [l, u]
    case (SOME((lower, true)), SOME((upper, true))) algorithm
      str := "Model error: Result of " + name + " outside the range "
        + realString(lower) + " <= " + ExpressionDump.printExpStr(rhs)
        + " <= " + realString(upper) + ". Unable to invert.";
      l := DAE.RELATION(DAE.RCONST(lower), DAE.LESSEQ(tp), rhs, -1, NONE());
      u:= DAE.RELATION(rhs, DAE.LESSEQ(tp), DAE.RCONST(upper), -1, NONE());
    then (str, DAE.LBINARY(l, DAE.AND(tp), u));

    // range [l, u)
    case (SOME((lower, true)), SOME((upper, false))) algorithm
      str := "Model error: Result of " + name + " outside the range "
        + realString(lower) + " <= " + ExpressionDump.printExpStr(rhs)
        + " < " + realString(upper) + ". Unable to invert.";
      l:= DAE.RELATION(DAE.RCONST(lower), DAE.LESSEQ(tp), rhs, -1, NONE());
      u:= DAE.RELATION(rhs, DAE.LESS(tp), DAE.RCONST(upper), -1, NONE());
    then (str, DAE.LBINARY(l, DAE.AND(tp), u));

    // range (l, u]
    case (SOME((lower, false)), SOME((upper, true))) algorithm
      str := "Model error: Result of " + name + " outside the range "
        + realString(lower) + " < " + ExpressionDump.printExpStr(rhs)
        + " <= " + realString(upper) + ". Unable to invert.";
      l:= DAE.RELATION(DAE.RCONST(lower), DAE.LESS(tp), rhs, -1, NONE());
      u:= DAE.RELATION(rhs, DAE.LESSEQ(tp), DAE.RCONST(upper), -1, NONE());
    then (str, DAE.LBINARY(l, DAE.AND(tp), u));

    // range (l, u)
    case (SOME((lower, false)), SOME((upper, false))) algorithm
      str := "Model error: Result of " + name + " outside the range "
        + realString(lower) + " < " + ExpressionDump.printExpStr(rhs)
        + " < " + realString(upper) + ". Unable to invert.";
      l:= DAE.RELATION(DAE.RCONST(lower), DAE.LESS(tp), rhs, -1, NONE());
      u:= DAE.RELATION(rhs, DAE.LESS(tp), DAE.RCONST(upper), -1, NONE());
    then (str, DAE.LBINARY(l, DAE.AND(tp), u));

    // range [l, inf)
    case (SOME((lower, true)), NONE()) algorithm
      str := "Model error: Result of " + name + " should be "
        + ExpressionDump.printExpStr(rhs) + " >= " + realString(lower)
        + ". Unable to invert.";
      l:= DAE.RELATION(DAE.RCONST(lower), DAE.LESSEQ(tp), rhs, -1, NONE());
    then (str, l);

    // range (l, inf)
    case (SOME((lower, true)), NONE()) algorithm
      str := "Model error: Result of " + name + " should be "
        + ExpressionDump.printExpStr(rhs) + " > " + realString(lower)
        + ". Unable to invert.";
      l:= DAE.RELATION(DAE.RCONST(lower), DAE.LESS(tp), rhs, -1, NONE());
    then (str, l);

    // range (-inf, u]
    case (NONE(), SOME((upper, true))) algorithm
      str := "Model error: Result of " + name + " should be "
        + ExpressionDump.printExpStr(rhs) + " <= " + realString(upper)
        + ". Unable to invert.";
      u:= DAE.RELATION(rhs, DAE.LESSEQ(tp), DAE.RCONST(upper), -1, NONE());
    then (str, u);

    // range (-inf, u)
    case (NONE(), SOME((upper, false))) algorithm
      str := "Model error: Result of " + name + " should be "
        + ExpressionDump.printExpStr(rhs) + " < " + realString(upper)
        + ". Unable to invert.";
      u:= DAE.RELATION(rhs, DAE.LESS(tp), DAE.RCONST(upper), -1, NONE());
    then (str, u);
  end match;

  algo := DAE.ALGORITHM_STMTS({DAE.STMT_ASSERT(cond, DAE.SCONST(msg), DAE.ASSERTIONLEVEL_ERROR, DAE.emptyElementSource)});
  assEq := BackendDAE.ALGORITHM(0, algo, DAE.emptyElementSource, DAE.EXPAND(), BackendDAE.EQ_ATTR_DEFAULT_UNKNOWN);
end makeDomainAssert;

protected function makeInitialGuess
  input DAE.Type tp;
  input DAE.Exp iExp1;
  input DAE.Exp iExp2;
  output DAE.Exp oExp;
protected
  DAE.Exp con, e;
algorithm
 con := Expression.makePureBuiltinCall("initial", {}, tp);
 e := Expression.traverseExpBottomUp(iExp2, makeInitialGuess2, (iExp1, "pre", tp, true));
 oExp := Expression.traverseExpBottomUp(iExp2, makeInitialGuess2, (iExp1, "pre", tp, false));
 oExp := DAE.IFEXP(con, e, oExp);
end makeInitialGuess;

protected function makeInitialGuess2
  input DAE.Exp iExp;
  input tuple<DAE.Exp, String, DAE.Type, Boolean> itpl;
  output DAE.Exp oExp;
  output tuple<DAE.Exp, String, DAE.Type, Boolean> otpl = itpl;
algorithm
  oExp := match(iExp, itpl)
    local
      DAE.ComponentRef cr1,cr2;
      DAE.Type tp;
      String fun;
      DAE.Exp e;

    case (DAE.CREF(componentRef=cr1), (DAE.CREF(componentRef=cr2), fun, tp, _))
      guard(ComponentReference.crefEqual(cr1, cr2)) algorithm
      e := Expression.makePureBuiltinCall(fun, {iExp}, tp);
    then e;

    case (_, (_, _, tp, true)) algorithm
      try
        SOME(e) := makeInitialGuess3(iExp, tp);
      else
        e := iExp;
      end try;
    then e;

    else iExp;
  end match;
end makeInitialGuess2;

protected function makeInitialGuess3
  input DAE.Exp iExp;
  input DAE.Type tp;
  output Option<DAE.Exp> oExp;
algorithm
  oExp := match(iExp)
          local DAE.Exp e, con, o;

          case(DAE.CALL(path = Absyn.IDENT(name = "log"), expLst={e}))
          equation
            con =  DAE.RELATION(e, DAE.LESSEQ(tp), DAE.RCONST(0.0), -1, NONE());
            o = DAE.IFEXP(con, DAE.RCONST(-1/0.000000001), iExp);
          then SOME(o);

          case(DAE.CALL(path = Absyn.IDENT(name = "log10"), expLst={e}))
          equation
            con =  DAE.RELATION(e, DAE.LESSEQ(tp), DAE.RCONST(0.0), -1, NONE());
            o = DAE.IFEXP(con, DAE.RCONST(-1/0.000000001), iExp);
          then SOME(o);

          case(DAE.CALL(path = Absyn.IDENT(name = "sqrt"), expLst={e}))
          equation
            con =  DAE.RELATION(e, DAE.LESSEQ(tp), DAE.RCONST(0.0), -1, NONE());
            o = DAE.IFEXP(con, DAE.RCONST(0.0), iExp);
          then SOME(o);

          case(DAE.BINARY(exp2=e))
          equation
            con =  DAE.RELATION(e, DAE.EQUAL(tp), DAE.RCONST(0.0), -1, NONE());
            o = DAE.IFEXP(con, DAE.RCONST(1.0), iExp);
          then SOME(o);

          else NONE();

         end match;

end makeInitialGuess3;

protected function helpInvCos
  input DAE.Exp acosy;
  input DAE.Exp x;
  input DAE.Type tp;
  input Boolean neg;
  output DAE.Exp k;
algorithm
  k := if neg then
         Expression.expAdd(x,acosy)
       else
         Expression.expSub(x,acosy);
  k := Expression.makeDiv(k, Expression.expMul(DAE.RCONST(2.0), DAE.PI));
  k := Expression.makePureBuiltinCall("$_round",{k},tp);

end helpInvCos;

protected function helpInvSin
  input DAE.Exp asiny;
  input DAE.Exp x;
  input DAE.Type tp;
  input Boolean neg;
  output DAE.Exp k;
algorithm
  k := if neg then
         Expression.expAdd(x,asiny)
       else
         Expression.expSub(x,asiny);
  k := Expression.makeDiv(k, Expression.expMul(DAE.RCONST(2.0), DAE.PI));
  if neg then
    k := Expression.expSub(k, DAE.RCONST(0.5));
  end if;
  k := Expression.makePureBuiltinCall("$_round",{k},tp);
end helpInvSin;

protected function helpInvCos2
  input DAE.Exp k;
  input DAE.Exp acosy;
  input DAE.Type tp;
  input Boolean neg;
  output DAE.Exp x;
algorithm

  x := if neg then Expression.negate(acosy) else acosy;
  x := Expression.expAdd(x, Expression.expMul(k, Expression.expMul(DAE.RCONST(2.0), DAE.PI)));

end helpInvCos2;

protected function helpInvSin2
  input DAE.Exp k;
  input DAE.Exp asiny;
  input DAE.Type tp;
  input Boolean neg;
  output DAE.Exp x;
protected
  DAE.Exp e;
algorithm

  x := if neg then Expression.negate(asiny) else asiny;
  e := Expression.expMul(k, Expression.expMul(DAE.RCONST(2.0), DAE.PI));
  e := if neg then Expression.expAdd(e, DAE.PI) else e;
  x := Expression.expAdd(x, e);

end helpInvSin2;

protected function helpInvCos3
  input DAE.Exp x1;
  input DAE.Exp x2;
  input DAE.Exp x;
  input DAE.Type tp;
  output DAE.Exp y;
protected
  DAE.Exp diffx1 = absDiff(x1,x,tp);
  DAE.Exp diffx2 = absDiff(x2,x,tp);
  DAE.Exp con = DAE.RELATION(diffx1, DAE.LESS(tp), diffx2, -1, NONE());
algorithm
  con := Expression.makeNoEvent(con);
  y := DAE.IFEXP(con, x1, x2);
end helpInvCos3;

protected function absDiff
  input DAE.Exp x;
  input DAE.Exp y;
  input DAE.Type tp;
  output DAE.Exp z;
algorithm
  z := Expression.expSub(x,y);
  z := Expression.makePureBuiltinCall("abs",{z},tp);
end absDiff;

annotation(__OpenModelica_Interface="backend");
end ExpressionSolve;
