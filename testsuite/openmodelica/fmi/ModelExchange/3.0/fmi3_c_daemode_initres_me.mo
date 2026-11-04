model DaeModeInitResC "Initial equations exported as residuals beside the fmi-ls-dae residuals"
  parameter Real k = 2;
  Real x1(start = 1) "initialized by the initial equations";
  Real x2(start = 1) "initialized by the initial equations";
  Real x3(start = 3, fixed = true) "a start value, not an initial residual";
  Real x4(start = 1) "steady state, through an auxiliary";
  Real y "in a linear loop with z";
  Real z(nominal = 0.1);
  Real w "an auxiliary variable of the DAE formulation";
equation
  der(x1) = y - x1 + 1;
  der(x2) = x1 - k*x2;
  der(x3) = -x3;
  der(x4) = w - x4;
  y + 2*z = x1;
  y - 3*z = x2;
  w = x1*z + x3;
initial equation
  x1 = 2*x2;
  der(x1) = 0;
  w = x4;
  annotation(
    __OpenModelica_commandLineOptions = "--daeMode",
    experiment(StopTime = 1, Tolerance = 1e-6));
end DaeModeInitResC;
