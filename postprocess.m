

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

%%
%---fit gaussians (this should really be in a test suite)
%-----we're only fitting rho_2 and rho_3
n_start_linear_fit = 20;
j = 1;
for i = 1:nt
    fit_ca_2 = fit(rr.',sol(i,:,2).','gauss1');
    fit_ca_3 = fit(rr.',sol(i,:,3).','gauss1');
    Gamma_t_2_a(j) = fit_ca_2.c1.^2/4;
%     Gamma_t_2_a_alt(j) = 1/fit_ca_2.a1.^2/(4*pi);
    Gamma_t_3_a(j) = fit_ca_3.c1.^2/4;
%     Gamma_t_3_a_alt(j) = (1/fit_ca_3.a1.^2)/(4*pi);
    j = j+1;
end
%---linear fit time
t_fit = t(n_start_linear_fit:end).';
%---fit linear poly to Gamma*t
fit_Gamma_2 = fit(t_fit,Gamma_t_2_a(n_start_linear_fit:end).','poly1');
fit_Gamma_3 = fit(t_fit,Gamma_t_3_a(n_start_linear_fit:end).','poly1');
% fit_Gamma_2_alt = fit(t_fit,Gamma_t_2_a_alt(n_start_linear_fit:end).','poly1');
% fit_Gamma_3_alt = fit(t_fit,Gamma_t_3_a_alt(n_start_linear_fit:end).','poly1');
%---plot 2
figure
plot(...
    t,Gamma_t_2_a,...
    'LineWidth',1 ...
); hold on
plot(...
    t,fit_Gamma_2(t),'--',...
    'LineWidth',2 ...
); hold on
% plot(...
%     t,Gamma_t_2_a_alt,...
%     'LineWidth',1 ...
% ); hold on
% plot(...
%     t,fit_Gamma_2_alt(t),'--',...
%     'LineWidth',2 ...
% );
grid on
xlabel('dimensionless time')
ylabel('nuclear distribution Gaussian variance')
legend(...
    'Gaussian variance','linear fit',...
    'Location','NorthWest'...
)
%---save
matlab2tikz('figures/gaussian_variance_vs_time_2.tex','width','\figurewidth')
%---plot 3
figure
plot(...
    t,Gamma_t_3_a,...
    'LineWidth',1 ...
); hold on
plot(...
    t,fit_Gamma_3(t),'--',...
    'LineWidth',2 ...
); hold on
% plot(...
%     t,Gamma_t_3_a_alt,...
%     'LineWidth',1 ...
% ); hold on
% plot(...
%     t,fit_Gamma_3_alt(t),'--',...
%     'LineWidth',2 ...
% );
grid on
xlabel('dimensionless time')
ylabel('electron distribution Gaussian variance')
legend(...
    'Gaussian variance','linear fit',...
    'Location','NorthWest'...
)
%---save
matlab2tikz('figures/gaussian_variance_vs_time_3.tex','width','\figurewidth')