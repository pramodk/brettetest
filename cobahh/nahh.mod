NEURON {
  SUFFIX nahh
  USEION na READ ena WRITE ina
  RANGE gmax
}

PARAMETER {
  gmax = 0.1 (mho/cm2)
}

ASSIGNED {
  v (millivolt)
  ena (millivolt)
  ina (milliamp/cm2)
  g (mho/cm2)
  am (/ms)
  bm (/ms)
  ah (/ms)
  bh (/ms)
}

STATE {
  m h
}

INITIAL {
  rates(v)
  m = am/(am + bm)
  h = ah/(ah + bh)
}

BREAKPOINT {
  SOLVE state METHOD cnexp
  g = gmax * m*m*m*h
  ina = g*(v - ena)
}

DERIVATIVE state {
  rates(v)
  m' = am*(1 - m) - bm*m
  h' = ah*(1 - h) - bh*h
}

PROCEDURE rates(v(millivolt)) {
  UNITSOFF
    am = 1.28*einstein(0.25*(v + 50))
    bm = 1.4*einstein(-0.2*(v + 23))
    ah = 0.128*exp(-0.05556*(v + 46))
    bh = 4/(1 + exp(0.2*(-23 - v)))
  UNITSON
}

FUNCTION einstein(x) {
  if (fabs(x) < 1e-6) {
    einstein =  1 + x/2
  }else{
    einstein = x/(1 - exp(-x))
  }
}
