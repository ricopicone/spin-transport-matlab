% --------------------------------------------------------------
function [cm,f,s] = sim_pde(r,t,p,DpDx)

% constants and parameters
global ii g G D c rho_3_langevin T1 mu_p mu_e B0 B_d kB temp T12 T13 freq duty rmax sim T

cm = [1; 1; 1];                                  
f = [...
    1+G; ...
    1; ...
    G ...
].*DpDx;

% T1 terms
switch sim
    case 'performed'
        simPulsePerformedExperiment;
    case 'rugar1'
        T12Term = 0;
        T13Term = 0;
    otherwise
        T12Term = 0;
        T13Term = 0;
end

% separation terms
s = [...
    -c^2/(1+D)*( (1-p(2)^2) + G*D*g^2*(1-p(3)^2) )*atanh(p(1))-c/(1+D)*(DpDx(2)-G*D*g*DpDx(3)); ...
    c*( (1-p(2)^2)/(1-p(1)^2)*DpDx(1) - 2*atanh(p(1))*p(2)*DpDx(2) ); ...
    -c*G*g*( (1-p(3)^2)/(1-p(1)^2)*DpDx(1) - 2*atanh(p(1))*p(3)*DpDx(3) ) ...
];  

ii = ii+1;
if mod(ii,10000)==0
    fprintf('%d  %g\n',ii,t)
end

end