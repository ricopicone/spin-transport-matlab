

%% Post-processing
rho_1 = sol(index_vec,:,1);  
rho_2 = sol(index_vec,:,2);  
rho_3 = sol(index_vec,:,3);
rho = sol(index_vec,:,:);
omega_1 = -(Delta_2+Delta_3)/Delta_2*atanh(sol(index_vec,:,1));  
omega_2 = -Delta_2/Delta_2*atanh(sol(index_vec,:,2));  
omega_3 = -Delta_3/Delta_2*atanh(sol(index_vec,:,3));
tdec = t(index_vec);

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
current = solcurrent(index_vec,:,:);

%---lambda (dynamic figure of merit)
kappa = -c.* ( 1 - sol(:,:,2).^2 ).* atanh( sol(:,:,1) ); % kappa really only defined in the steady state
[lambdaFull, dtkappaFull] = gradient(kappa, dr, dt);
lambda = lambdaFull(index_vec,:,:);

%---lambda cumulative
lambdaFullCum = dt*cumtrapz(lambdaFull);
lambdaCum = lambdaFullCum(index_vec,:);