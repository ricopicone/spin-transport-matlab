function out = B1(t,freq,duty,T_s_ramp,B1max)

% turn off B1 for a species by setting B1max to 0

td = mod(t,1/freq);


if td < T_s_ramp/2
  out = B1max*s_ramp_up(t,1/T_s_ramp);
elseif td < duty/freq
  out = B1max;
elseif td < duty/freq + T_s_ramp/2
  out = B1max*s_ramp_dn(t,1/T_s_ramp,duty,1/freq);
else
  out = 0;
end

end