model DaeModeInputDerC "Index 2: a state tied to an input"
  input Real u(start = 1);
  input Real v(start = 0) "its derivative is not needed";
  Real x(start = 1) "x = u: der(x) = der(u)";
  Real y "determined by der(x), i.e. by der(u)";
  Real s(start = 0, fixed = true);
  Real a(start = 1) "differentiated by index reduction";
  Real b "determined by der(a)";
equation
  x = u;
  der(x) = y - x;
  der(s) = y - s + v;
  a^2 = u + 1;
  der(a) = b;
  annotation(
    __OpenModelica_commandLineOptions = "--daeMode",
    experiment(StopTime = 1, Tolerance = 1e-6));
end DaeModeInputDerC;

model DaeModeInputDerWrittenC "der(u) written in the model"
  extends DaeModeInputDerC;
  Real w = 2*der(u);
end DaeModeInputDerWrittenC;
