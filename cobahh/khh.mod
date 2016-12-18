NEURON {
  SUFFIX khh
  USEION k READ ek WRITE ik
  RANGE gmax
}

PARAMETER {
  gmax = 0.03 (mho/cm2)
}

ASSIGNED {
  v (millivolt)
  ek (millivolt)
  ik (milliamp/cm2)
  g (mho/cm2)
  an (/ms)
  bn (/ms)
}

STATE {
  n
}

INITIAL {
  rates(v)
  n = an/(an + bn)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  g = gmax * n*n*n*n
  ik = g*(v - ek)
}

DERIVATIVE state {
  rates(v)
  n' = an*(1 - n) - bn*n
}

PROCEDURE rates(v(millivolt)) {
  UNITSOFF
    an = 0.16*einstein(0.2*(v + 48))
    bn = 0.5*exp(-0.025*(v + 53))
  UNITSON
}

FUNCTION einstein(x) {
  if (fabs(x) < 1e-6) {
    einstein =  1 + x/2
  }else{
    einstein = x/(1 - exp(-x))
  }
}
