%% Spatial grid
if mod(nr,2)  == 0
    nr = nr + 1;
end
dr=2*rmax/(nr-1);
for i=1:nr
    rr(i)=(i-1)*dr - rmax;
end

u0=PiconeTwoSpeciesTransport_ic(rmax);

%% Normalized temporal grid
%decim=10;  % decimate the time solutions
m = 0;
dt = T/(nt-1);
t=0:dt:T;