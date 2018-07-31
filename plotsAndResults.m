%% Results
p2slope_t=mu_p*B_d/kB/temp;     %---theoretical slope
fprintf('p2slope_t             = %g\n',p2slope_t);
p2slope_s=(rho_2(end,end)-rho_2(end,1))/(rr(end)-rr(1)); %---simlation slope 
fprintf('p2slope_s             = %g\n',p2slope_s);
p3slope_t=mu_e*B_d/kB/temp;     %---theoretical slope 
fprintf('p3slope_t             = %g\n',p3slope_t);
p3slope_s=(rho_3(end,end)-rho_3(end,1))/(rr(end)-rr(1)); %---simlation slope
fprintf('p3slope_s             = %g\n',p3slope_s);

%% Plots: 3D
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4) scrsz(3) 2/3*scrsz(4)])
subplot(131)
surf(tdec,rr,rho_1.','FaceColor','blue','FaceAlpha',0.01,'EdgeColor','interp')
axis vis3d
subplot(132)
surf(tdec,rr,rho_2.','FaceColor','blue','FaceAlpha',0.01,'EdgeColor','interp')
axis vis3d
subplot(133)
surf(tdec,rr,rho_3.','FaceColor','blue','FaceAlpha',0.01,'EdgeColor','interp')
axis vis3d

%% Plots: 2D


grey = [ 0.5, 0.5, 0.5 ];
C = linspecer(nTraces); 
% C = zeros( nTraces, 3 );
% for i = 1:nTraces % blue to red
%     C(i,1) = ( i - 1 ) / ( nTraces - 1 );
%     C(i,3) = ( nTraces - i ) / ( nTraces - 1 );
% end

lwidth = 1;
%---rho
figure('Position',[1 scrsz(4) scrsz(3) 2/3*scrsz(4)])
for i = 1:3
    subplot(2,3,i)
    for j = 1:nTraces
        switch sim
            case 'rugar1'
                plot([rr(1),rr(end)],[0,0],'LineWidth',0.5,'Color',grey); hold on
        end
        plot(rr, rho(j,:,i), 'LineWidth', lwidth,'Color',C(j,:)); hold on
        xlim([-range,range]+edge+range);
        title(['$$\rho_',num2str(i),'$$'],'Interpreter','Latex')
    end
    myLims(1,i,:) = ylim;
end
for i = 1:3
    subplot(2,3,(i+3))
    for j = 1:nTraces
        switch sim
            case 'rugar1'
                plot([rr(1),rr(end)],[0,0],'LineWidth',0.5,'Color',grey); hold on
        end
        plot(rr, current(j,:,i), 'LineWidth', lwidth,'Color',C(j,:)); hold on
        xlim([-range,range]+edge+range);
        title(['$$j_',num2str(i),'$$'],'Interpreter','Latex')
        xlabel('$$\overline{r} = \frac{\partial_r B|_{r = 0}}{B_d} r$$','Interpreter','Latex')
    end
    myLims(2,i,:) = ylim;
end


% figure;
% colormap(C);
% cbar = colorbar('North','XTickLabel',{'t = 0','t = 1'},'XMinorTick','off','XTick',[1,length(C)+1]);

if saveFlag == 1
    switch sim
        case 'rugar1'
            matlab2tikz('figures/proposedTikz.tex','width','\figurewidth')
        case 'performed'
            matlab2tikz('figures/performedTikz.tex','width','\figurewidth')
    end
end

%% lambda
figure
subplot(2,1,1)
for j = 1:nTraces
    plot(rr, lambda(j,:), 'LineWidth', lwidth,'Color',C(j,:)); hold on
    xlim([-range,range]+edge+range);
    title('$$\lambda(\overline{t},\overline{r})$$','Interpreter','Latex')
    %xlabel('$$\overline{r} = \frac{\partial_r B|_{r = 0}}{B_d} r$$','Interpreter','Latex')
end
subplot(2,1,2)
for j = 1:nTraces
    plot(rr, lambdaCum(j,:), 'LineWidth', lwidth,'Color',C(j,:)); hold on
    xlim([-range,range]+edge+range);
    title('$$\int_0^{\overline{t}}\lambda(\overline{t}^\prime,\overline{r}) d\overline{t}^\prime$$','Interpreter','Latex')
    xlabel('$$\overline{r} = \frac{\partial_r B|_{r = 0}}{B_d} r$$','Interpreter','Latex')
end

if saveFlag == 1
    switch sim
        case 'rugar1'
            matlab2tikz('figures/proposedlambdaTikz.tex','width','\figurewidth','height','\figureheight')
        case 'performed'
            matlab2tikz('figures/performedlambdaTikz.tex','width','\figurewidth','height','\figureheight')
    end
end

%% kappa
% figure
% plot(rr, kappa(1,:), 'LineWidth', lwidth,'Color',C(1,:)); hold on
% xlim([-range,range]+edge+range);
% title('$$\kappa$$','Interpreter','Latex')
% xlabel('$$\overline{r} = \frac{\partial_r B|_{r = 0}}{B_d} r$$','Interpreter','Latex')

%% Plots: 2D Animation

saveAnimation = 0;
nFrames = 60;
grey = [ 0.5, 0.5, 0.5 ]; 
indexVec2 = round(logspace(0,log10(length(t)),nFrames));
C = linspecer(nFrames);

lwidth = 1;
%---rho
for j = 1:nFrames
    close
    figure('Position',[1 scrsz(4) scrsz(3) 2/3*scrsz(4)])
    for i = 1:3
        subplot(2,3,i)
            switch sim
                case 'rugar1'
                    plot([rr(1),rr(end)],[0,0],'LineWidth',0.5,'Color',grey); hold on
            end
            plot(rr, sol(indexVec2(j),:,i), 'LineWidth', lwidth,'Color',C(j,:)); hold on
            xlim([-range,range]+edge+range);
            title(['$$\rho_',num2str(i),'$$'],'Interpreter','Latex')
            xlabel('$$\overline{r}$$','Interpreter','Latex')
            ylim(myLims(1,i,:))
    end
    for i = 1:3
        subplot(2,3,(i+3))
            switch sim
                case 'rugar1'
                    plot([rr(1),rr(end)],[0,0],'LineWidth',0.5,'Color',grey); hold on
            end
            plot(rr, solcurrent(indexVec2(j),:,i), 'LineWidth', lwidth,'Color',C(j,:)); hold on
            xlim([-range,range]+edge+range);
            title(['$$j_',num2str(i),'$$'],'Interpreter','Latex')
%             xlabel('$$\overline{r} = \frac{\partial_r B|_{r = 0}}{B_d} r$$','Interpreter','Latex')
            xlabel('$$\overline{r}$$','Interpreter','Latex')
            ylim(myLims(2,i,:))
    end

    if ( ( saveFlag == 1 ) && ( saveAnimation == 1 ) )
        switch sim
            case 'rugar1'
                matlab2tikz(['figures/proposedAnimationTikz',num2str(j),'.tex'],'width','\figurewidth')
            case 'performed'
                matlab2tikz(['figures/performedAnimationTikz',num2str(j),'.tex'],'width','\figurewidth')
        end
    end
end