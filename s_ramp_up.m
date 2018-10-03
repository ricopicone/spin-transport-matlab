function out = s_ramp_up(td,freq)
  out = (1-cos(2*pi*freq*td))/2;