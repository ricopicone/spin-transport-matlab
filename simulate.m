indexVec = floor(logspace(0,log10(length(t)),nTraces));
ii=0;
options = ...
    odeset( ...
        'AbsTol',1e-8, ...
        'RelTol',1e-5 ...
    );
sol =   pdepe(m, ...
        @PiconeTwoSpeciesTransport_pde, ...
        @PiconeTwoSpeciesTransport_ic, ...
        @PiconeTwoSpeciesTransport_bc, ...
        rr, ...
        t, ...
        options);
    
%% Post-processing
rho_1 = sol(indexVec,:,1);  
rho_2 = sol(indexVec,:,2);  
rho_3 = sol(indexVec,:,3);
rho = sol(indexVec,:,:);
omega_1 = -(Delta_2+Delta_3)/Delta_2*atanh(sol(indexVec,:,1));  
omega_2 = -Delta_2/Delta_2*atanh(sol(indexVec,:,2));  
omega_3 = -Delta_3/Delta_2*atanh(sol(indexVec,:,3));
tdec = t(indexVec);

%---current
dtsol = zeros(size(sol));
drsol = zeros(size(sol));
for i = 1:3
    [drsol(:,:,i), dtsol(:,:,i)] = gradient(sol(:,:,i), dr, dt);
end

solcurrent = zeros(size(sol));
solcurrent(:,:,1) = -(1+G).*drsol(:,:,1);
solcurrent(:,:,2) = -c.*( 1 - sol(:,:,2).^2 ).* atanh(sol(:,:,1)) - drsol(:,:,2);
solcurrent(:,:,3) = G.* ( g*c.*( 1 - sol(:,:,3).^2 ).* atanh(sol(:,:,1)) - drsol(:,:,3) );
current = solcurrent(indexVec,:,:);

%---lambda (dynamic figure of merit)
kappa = -c.* ( 1 - sol(:,:,2).^2 ).* atanh( sol(:,:,1) ); % kappa really only defined in the steady state
[lambdaFull, dtkappaFull] = gradient(kappa, dr, dt);
lambda = lambdaFull(indexVec,:,:);

%---lambda cumulative
lambdaFullCum = dt*cumtrapz(lambdaFull);
lambdaCum = lambdaFullCum(indexVec,:);