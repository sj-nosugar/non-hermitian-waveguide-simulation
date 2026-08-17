clear;clc;tic;
%% prm
delta = 0.5;
lam0=1.55E-6; 
wco=0.25*lam0;
%% csl.data
csl.da35=importdata("nco35-1.txt");%nco=3.5,result,rco=1[um]  
% csl.da35dif=importdata("rsl-txt/nco35005_eff_difrco.txt");%nco=3.5005,result with different rco
wg.distance=str2double(csl.da35.textdata(:,1));%
wg.beta=csl.da35.data(:,1);

%% prepare
wg.kappa=zeros(length(wg.beta)/2,1);
wg.dis=zeros(length(wg.beta)/2,1);


for ii =1:length(wg.distance)/2
    wg.kappa(ii)=sqrt((wg.beta(2*ii)-wg.beta(2*ii-1))^2/4);
    wg.dis(ii)=wg.distance(2*ii);
end

%% fit the result and compare
wg.t=(min(wg.kappa)+max(wg.kappa))/2;
% lsqcurevefit
fitf1 = @(p,x) p(1)*exp(p(2)*x)+p(3);
fit1.p0=[620000,-13000000,85];
fit1.p=lsqcurvefit(fitf1,fit1.p0,wg.dis,wg.kappa);
fit1.d=linspace(min(wg.dis),max(wg.dis),1000);
fit1.k=fitf1(fit1.p,fit1.d);
% % ln(wg.kappa)=ln(a)+b*x
% wg.kappa_lg=log(wg.kappa);
% fit2.p=polyfit(wg.dis,wg.kappa_lg,1);
% fit2.aln=exp(fit2.p(2));
% fit2.bln=fit2.p(1);
% fit2.k=fit2.aln*exp(fit2.bln*fit1.d);

wg.kappa=wg.kappa/wg.t;

% origin result and curve result.
figure
plot(wg.dis/wco, wg.kappa, 'ro', 'MarkerSize', 8, 'DisplayName', 'origin');
hold on;
plot(fit1.d/wco, fit1.k/wg.t, 'b-', 'LineWidth', 2, 'DisplayName', 'curvefit');
% plot(fit1.d,fit2.k/wg.t,'g--','LineWidth',2,'DisplayName','polyfit')
grid on;
legend('Location', 'best');
title('指数函数拟合结果','FontSize',16);
xlabel('distance(x[m])/wco[m]','FontSize',14);
ylabel('coupling constant \kappa','FontSize',14);
axis tight;
% % 绘制残差图
% subplot(2,1,2);
% cp.kappa_pred = fitf1(fit1.p, wg.dis);
% cp.residuals = wg.kappa - cp.kappa_pred;
% stem(wg.dis, cp.residuals, 'filled', 'MarkerSize', 6);
% grid on;
% title('残差图');
% xlabel('x');
% ylabel('残差');
% line([min(wg.dis), max(wg.dis)], [0, 0], 'Color', 'r', 'LineStyle', '--');
% cp.SSR = sum((cp.kappa_pred - mean(wg.kappa)).^2);
% cp.SST = sum((wg.kappa - mean(wg.kappa)).^2);
% cp.R_squared = cp.SSR / cp.SST;

% 输出拟合参数和R²
fprintf('拟合参数: a=%.4f, b=%.4f, c=%.4f\n', fit1.p(1), fit1.p(2), fit1.p(3));
fprintf('耦合常数量纲值 = %.4f\n',wg.t);
% fprintf('拟合优度R²: %.4f\n', cp.R_squared);
%% choose the t and $\delta$
coupling_function=@(x) fit1.p(1)*exp(fit1.p(2)*x)+fit1.p(3);
wg.tplus=(1+delta)*wg.t;
wg.tminus=(1-delta)*wg.t;

dis.xplus=fzero(@(x)coupling_function(x)-wg.tplus,wco)/wco;
dis.x0=fzero(@(x)coupling_function(x)-wg.t,wco)/wco;
dis.xminus=fzero(@(x)coupling_function(x)-wg.tminus,wco)/wco;
fprintf('dis_plus =%.4f, dis_minus = %.4f, dis_neutral=%.4f\n',...
    dis.xplus, dis.xminus, dis.x0);
fprintf('dis_plus[um] =%.4f, dis_minus[um] = %.4f, dis_neutral[um]=%.4f\n',...
    dis.xplus*wco*1e6, dis.xminus*wco*1e6, dis.x0*wco*1e6);

% wg.t=(max(wg.kappa)+min(wg.kappa))/2;
% wg.delta=delta;
% wg.tplusd=wg.t+wg.delta;
% wg.tminusd=wg.t-wg.delta;
% wg.dis0=mean(wg.dis);
% options=optimset('Display','off');
% rsl.dis_t=fzero(@(x) f(x)-wg.t,wg.dis0,options);
% t=rsl.dis_t/rco;
% rsl.dis_tplusd=fzero(@(x) f(x)-wg.tplusd,wg.dis0,options);
% t_pd=rsl.dis_tplusd/rco;
% rsl.dis_tminusd=fzero(@(x) f(x)-wg.tminusd,wg.dis0,options);
% t_md=rsl.dis_tminusd/rco;
% %% show the result
% figure('Position', [100, 100, 800, 600]);
% 
% % 绘制拟合曲线和目标点
% plot(fit1.d, fit1.k, 'b-', 'LineWidth', 2, 'DisplayName', 'fit result');
% hold on;
% 
% % 绘制水平参考线
% line([min(fit1.d), max(fit1.d)], [wg.t, wg.t], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
% line([min(fit1.d), max(fit1.d)], [wg.tplusd, wg.tplusd], 'Color', 'g', 'LineStyle', '--', 'LineWidth', 1);
% line([min(fit1.d), max(fit1.d)], [wg.tminusd, wg.tminusd], 'Color', 'm', 'LineStyle', '--', 'LineWidth', 1);
% 
% % 绘制目标点
% plot(rsl.dis_t, wg.t, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName',...
%     ['t对应点 (', num2str(rsl.dis_t, '%.4f'), ', ', num2str(wg.t, '%.4f'), ')']);
% plot(rsl.dis_tplusd, wg.tplusd, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName',...
%     ['t+δ对应点 (', num2str(rsl.dis_tplusd, '%.4f'), ', ', num2str(wg.tplusd, '%.4f'), ')']);
% plot(rsl.dis_tminusd, wg.tminusd, 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName',...
%     ['t-δ对应点 (', num2str(rsl.dis_tminusd, '%.4f'), ', ', num2str(wg.tminusd, '%.4f'), ')']);
% 
% grid on;
% legend('Location', 'best');
% title('指数函数特定值对应点');
% xlabel('distance x/[m]');
% ylabel('propagation constant \kappa');
% 
% format long e
% fprintf('delta=%.4f, x_{t}=%.4f, x_{t+delta}=%.4f, x_{t-delta}=%.4f', delta, rsl.dis_t, rsl.dis_tplusd, rsl.dis_tminusd);
% %%
%%
toc
