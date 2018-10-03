function out = s_ramp_dn(td,freq,duty,T)
  out = (1+cos(2*pi*freq*(td - duty*T)))/2;