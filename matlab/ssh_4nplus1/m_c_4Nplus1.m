% main.m - 主脚本，演示模型c的完整分析流程

%% 参数设置
clear; clc;close all;tic;
N = 20;          % 原胞数
gamma = 2.5;    % 增益损耗强度
J1 = 0.5;       % 原胞内跃迁强度
J2 = 1.5;       % 原胞间跃迁强度
t_max = 100;     % 演化时间
dt = 0.01;       % 时间步长
% printf setting
fontName = 'Helvetica';
fontSize = 16;
lineWidth = 3;

%% 1. 构建哈密顿量
fprintf('构建哈密顿量 (N=%d)...\n', N);
H = buildHamiltonian_c(N, gamma, J1, J2);
total_sites = 4*N + 1;  % 更新：4N+1个格点
if N <= 20
    fprintf('哈密顿量矩阵 (前10×10):\n');
    disp(full(H(1:min(10,total_sites), 1:min(10,total_sites))));
end

%% 2. 求解本征能谱
fprintf('求解本征值...\n');
[eigenvalues, eigenvectors] = solveSpectrum_c(H, N);

% 绘制能谱图
figure('Name', '本征值能谱', 'Position', [100 100 800 600]);
set(gca, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', lineWidth);
plot(real(eigenvalues), 'b.', 'MarkerSize', 15);
xlabel('本征态索引'); ylabel('能量实部');
title(sprintf('本征值实部 (N=%d, 总格点数=%d)', N, total_sites));
grid on; box on; axis tight;
print('-dpng','-r600',sprintf('Modelc_spectrum_N%d_gamma%.2f_J1%.2f_J2%.2f.png', N, gamma, J1, J2));


% 绘制复平面上的本征值分布
figure('Name', '复平面本征值分布', 'Position', [100 100 600 600]);
set(gca, 'FontName', fontName, 'FontSize', fontSize, 'LineWidth', lineWidth);
plot(real(eigenvalues), imag(eigenvalues), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
xlabel('实部'); ylabel('虚部');
title(sprintf('复平面上的本征值分布 (N=%d, 总格点数=%d)', N, total_sites),'FontName',fontName,'FontSize',fontSize);
grid on; box on;
axis tight;
print('-dpng','-r600',sprintf('Modelc_real_imag_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));

% %% 3. 设置三种初态
% fprintf('设置三种初态...\n');
% 
% % 计算特殊格点索引
% b_N_index = 2*N+1;      % b_N的索引
% a_Nplus1_index = 2*N+2;  % a_{N+1}的索引
% 
% % 初始化三种初态
% psi0_1 = zeros(total_sites, 1);  % 初态1: b_N上为1
% psi0_1(b_N_index) = 1;
% 
% psi0_2 = zeros(total_sites, 1);  % 初态2: a_{N+1}上为1
% psi0_2(a_Nplus1_index) = 1;
% 
% psi0_3 = zeros(total_sites, 1);  % 初态3: b_N和a_{N+1}上各为0.5
% psi0_3(b_N_index) = 0.5;
% psi0_3(a_Nplus1_index) = 0.5;
% psi0_3 = psi0_3 / norm(psi0_3);  % 归一化
% 
% %% 4. 时间演化
% fprintf('执行时间演化...\n');
% 
% % 模拟三种初态的时间演化
% [t1, psi_t1, prob_b_N1, ipr_evolution1] = timeEvolve_c(H, psi0_1, t_max, dt);
% [t2, psi_t2, prob_b_N2, ipr_evolution2] = timeEvolve_c(H, psi0_2, t_max, dt);
% [t3, psi_t3, prob_b_N3, ipr_evolution3] = timeEvolve_c(H, psi0_3, t_max, dt);
% 
% 
% 
% % 计算概率密度矩阵
% prob_density1 = abs(psi_t1).^2;
% prob_density2 = abs(psi_t2).^2;
% prob_density3 = abs(psi_t3).^2;
% 
% %% 5. 可视化演化过程
% 
% % 创建三个独立的演化热图
% 
% % 初态1的热图
% figure('Name', '初态1: b_N上为1 - 概率密度演化热图', 'Position', [100 100 1200 500]);
% imagesc(1:total_sites, t1, prob_density1');
% colorbar;
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% set(colorbar, 'FontName', fontName, 'FontSize', fontSize-2);
% ylabel('时间 t', 'FontName', fontName, 'FontSize', fontSize);
% xlabel('格点位置', 'FontName', fontName, 'FontSize', fontSize);
% title(sprintf('初态1: b_N上为1 (N=%d, 总格点数=%d)', N, total_sites), 'FontName', fontName, 'FontSize', fontSize);
% axis xy;
% grid on; box on;
% 
% % 标记特殊格点
% hold on;
% plot(b_N_index, 0, 'w.', 'MarkerSize', 15);
% text(b_N_index, -0.5, 'b_N', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', fontName, 'FontSize', fontSize);
% plot(a_Nplus1_index, 0, 'w.', 'MarkerSize', 15);
% text(a_Nplus1_index, -0.5, 'a_{N+1}', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', fontName, 'FontSize', fontSize);
% hold off;
% print('-dpng','-r600',sprintf('Modelc_density_initial1_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));
% 
% % 初态2的热图
% figure('Name', '初态2: a_{N+1}上为1 - 概率密度演化热图', 'Position', [100 100 1200 500]);
% imagesc(1:total_sites, t2, prob_density2');
% colorbar;
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% set(colorbar, 'FontName', fontName, 'FontSize', fontSize-2);
% ylabel('时间 t', 'FontName', fontName, 'FontSize', fontSize);
% xlabel('格点位置', 'FontName', fontName, 'FontSize', fontSize);
% title(sprintf('初态2: a_{N+1}上为1 (N=%d, 总格点数=%d)', N, total_sites), 'FontName', fontName, 'FontSize', fontSize+2);
% axis xy;
% grid on; box on;
% 
% % 标记特殊格点
% hold on;
% plot(b_N_index, 0, 'w.', 'MarkerSize', 15);
% text(b_N_index, -0.5, 'b_N', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', fontName, 'FontSize', fontSize);
% plot(a_Nplus1_index, 0, 'w.', 'MarkerSize', 15);
% text(a_Nplus1_index, -0.5, 'a_{N+1}', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', fontName, 'FontSize', fontSize);
% hold off;
% print('-dpng','-r600',sprintf('Modelc_density_initial2_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));
% 
% % 初态3的热图
% figure('Name', '初态3: b_N和a_{N+1}各0.5 - 概率密度演化热图', 'Position', [100 100 1200 500]);
% imagesc(1:total_sites, t3, prob_density3');
% colorbar;
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% set(colorbar, 'FontName', fontName, 'FontSize', fontSize-2);
% ylabel('时间 t', 'FontName', fontName, 'FontSize', fontSize);
% xlabel('格点位置', 'FontName', fontName, 'FontSize', fontSize);
% title(sprintf('初态3: b_N和a_{N+1}各0.5 (N=%d, 总格点数=%d)', N, total_sites), 'FontName', fontName, 'FontSize', fontSize+2);
% axis xy;
% grid on; box on;
% 
% % 标记特殊格点
% hold on;
% plot(b_N_index, 0, 'w.', 'MarkerSize', 15);
% text(b_N_index, -0.5, 'b_N', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', fontName, 'FontSize', fontSize);
% plot(a_Nplus1_index, 0, 'w.', 'MarkerSize', 15);
% text(a_Nplus1_index, -0.5, 'a_{N+1}', 'Color', 'w', 'HorizontalAlignment', 'center', 'FontName', fontName, 'FontSize', fontSize);
% hold off;
% print('-dpng','-r600',sprintf('Modelc_density_initial3_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));
% 
% %% 6. 最终概率分布对比
% figure('Name', '三种初态的最终概率分布对比', 'Position', [100 100 1400 800]);
% 
% % 初态1的最终分布
% subplot(3, 1, 1);
% bar(1:total_sites, prob_density1(:, end), 'FaceColor', [0.2, 0.6, 0.8]);
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% xlabel('格点位置', 'FontName', fontName, 'FontSize', fontSize);
% ylabel('概率', 'FontName', fontName, 'FontSize', fontSize);
% title('初态1: b_N上为1', 'FontName', fontName, 'FontSize', fontSize+2);
% grid on; box on;
% xlim([0, total_sites + 1]);
% ylim([0, max(prob_density1(:, end))*1.1]);
% hold on;
% plot(b_N_index, prob_density1(b_N_index, end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% plot(a_Nplus1_index, prob_density1(a_Nplus1_index, end), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
% text(b_N_index, prob_density1(b_N_index, end), '  b_N', 'Color', 'r', 'FontName', fontName, 'FontSize', fontSize);
% text(a_Nplus1_index, prob_density1(a_Nplus1_index, end), '  a_{N+1}', 'Color', 'g', 'FontName', fontName, 'FontSize', fontSize);
% hold off;
% 
% % 初态2的最终分布
% subplot(3, 1, 2);
% bar(1:total_sites, prob_density2(:, end), 'FaceColor', [0.8, 0.4, 0.2]);
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% xlabel('格点位置', 'FontName', fontName, 'FontSize', fontSize);
% ylabel('概率', 'FontName', fontName, 'FontSize', fontSize);
% title('初态2: a_{N+1}上为1', 'FontName', fontName, 'FontSize', fontSize+2);
% grid on; box on;
% xlim([0, total_sites + 1]);
% ylim([0, max(prob_density2(:, end))*1.1]);
% hold on;
% plot(b_N_index, prob_density2(b_N_index, end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% plot(a_Nplus1_index, prob_density2(a_Nplus1_index, end), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
% text(b_N_index, prob_density2(b_N_index, end), '  b_N', 'Color', 'r', 'FontName', fontName, 'FontSize', fontSize);
% text(a_Nplus1_index, prob_density2(a_Nplus1_index, end), '  a_{N+1}', 'Color', 'g', 'FontName', fontName, 'FontSize', fontSize);
% hold off;
% 
% % 初态3的最终分布
% subplot(3, 1, 3);
% bar(1:total_sites, prob_density3(:, end), 'FaceColor', [0.5, 0.3, 0.7]);
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% xlabel('格点位置', 'FontName', fontName, 'FontSize', fontSize);
% ylabel('概率', 'FontName', fontName, 'FontSize', fontSize);
% title('初态3: b_N和a_{N+1}各0.5', 'FontName', fontName, 'FontSize', fontSize+2);
% grid on; box on;
% xlim([0, total_sites + 1]);
% ylim([0, max(prob_density3(:, end))*1.1]);
% hold on;
% plot(b_N_index, prob_density3(b_N_index, end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% plot(a_Nplus1_index, prob_density3(a_Nplus1_index, end), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
% text(b_N_index, prob_density3(b_N_index, end), '  b_N', 'Color', 'r', 'FontName', fontName, 'FontSize', fontSize);
% text(a_Nplus1_index, prob_density3(a_Nplus1_index, end), '  a_{N+1}', 'Color', 'g', 'FontName', fontName, 'FontSize', fontSize);
% hold off;
% 
% sgtitle(sprintf('三种初态的最终概率分布对比 (N=%d, t=%.1f)', N, t_max), 'FontSize', 16, 'FontWeight', 'bold');
% print('-dpng','-r600',sprintf('Modelc_final_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));
% 
% %% 7. 目标格点概率演化对比
% figure('Name', '目标格点概率演化对比', 'Position', [100 100 1200 800]);
% 
% % 提取目标格点概率
% prob_b_N_1 = prob_density1(b_N_index, :);
% prob_a_Nplus1_1 = prob_density1(a_Nplus1_index, :);
% 
% prob_b_N_2 = prob_density2(b_N_index, :);
% prob_a_Nplus1_2 = prob_density2(a_Nplus1_index, :);
% 
% prob_b_N_3 = prob_density3(b_N_index, :);
% prob_a_Nplus1_3 = prob_density3(a_Nplus1_index, :);
% 
% % b_N格点概率演化
% subplot(2, 1, 1);
% plot(t1, prob_b_N_1, 'LineWidth', 2, 'Color', [0.2, 0.6, 0.8], 'DisplayName', '初态1');
% hold on;
% plot(t2, prob_b_N_2, 'LineWidth', 2, 'Color', [0.8, 0.4, 0.2], 'DisplayName', '初态2');
% plot(t3, prob_b_N_3, 'LineWidth', 2, 'Color', [0.5, 0.3, 0.7], 'DisplayName', '初态3');
% hold off;
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% xlabel('时间 t', 'FontName', fontName, 'FontSize', fontSize);
% ylabel('概率 |ψ(b_N)|²', 'FontName', fontName, 'FontSize', fontSize);
% title('b_N格点上的概率演化', 'FontName', fontName, 'FontSize', fontSize+2);
% grid on; box on;
% legend('Location', 'best', 'FontName', fontName, 'FontSize', fontSize);
% 
% % a_{N+1}格点概率演化
% subplot(2, 1, 2);
% plot(t1, prob_a_Nplus1_1, 'LineWidth', 2, 'Color', [0.2, 0.6, 0.8], 'DisplayName', '初态1');
% hold on;
% plot(t2, prob_a_Nplus1_2, 'LineWidth', 2, 'Color', [0.8, 0.4, 0.2], 'DisplayName', '初态2');
% plot(t3, prob_a_Nplus1_3, 'LineWidth', 2, 'Color', [0.5, 0.3, 0.7], 'DisplayName', '初态3');
% hold off;
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% xlabel('时间 t', 'FontName', fontName, 'FontSize', fontSize);
% ylabel('概率 |ψ(a_{N+1})|²', 'FontName', fontName, 'FontSize', fontSize);
% title('a_{N+1}格点上的概率演化', 'FontName', fontName, 'FontSize', fontSize+2);
% grid on; box on;
% legend('Location', 'best', 'FontName', fontName, 'FontSize', fontSize);
% print('-dpng','-r600',sprintf('Modelc_prob_compare_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));
% 
% 
% 
% %% 8. IPR演化对比
% figure('Name', 'IPR演化对比', 'Position', [100 100 1200 600]);
% 
% % 绘制三种初态的IPR演化
% subplot(1, 2, 1);
% plot(t1, ipr_evolution1, 'LineWidth', 2, 'Color', [0.2, 0.6, 0.8], 'DisplayName', '初态1: b_N上为1');
% hold on;
% plot(t2, ipr_evolution2, 'LineWidth', 2, 'Color', [0.8, 0.4, 0.2], 'DisplayName', '初态2: a_{N+1}上为1');
% plot(t3, ipr_evolution3, 'LineWidth', 2, 'Color', [0.5, 0.3, 0.7], 'DisplayName', '初态3: b_N和a_{N+1}各0.5');
% hold off;
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% xlabel('时间 t', 'FontName', fontName, 'FontSize', fontSize);
% ylabel('IPR', 'FontName', fontName, 'FontSize', fontSize);
% title('三种初态的IPR演化对比', 'FontName', fontName, 'FontSize', fontSize+2);
% grid on; box on;
% legend('Location', 'best', 'FontName', fontName, 'FontSize', fontSize);
% xlim([0, t_max]);
% ylim([0, 1]);
% 
% % 添加理论极限线
% hold on;
% plot([0, t_max], [1/total_sites, 1/total_sites], 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('均匀极限 1/N=%.4f', 1/total_sites));
% plot([0, t_max], [1, 1], 'r--', 'LineWidth', 1.5, 'DisplayName', '完全局域极限');
% hold off;
% 
% % 绘制最终IPR值的柱状图对比
% subplot(1, 2, 2);
% final_ipr = [ipr_evolution1(end), ipr_evolution2(end), ipr_evolution3(end)];
% bar_colors = {[0.2, 0.6, 0.8], [0.8, 0.4, 0.2], [0.5, 0.3, 0.7]};
% bar(1:3, final_ipr, 'FaceColor', 'flat');
% for i = 1:3
%     bar_handles = findobj(gca, 'Type', 'bar');
%     bar_handles.CData(i, :) = bar_colors{i};
% end
% 
% % 添加数值标签
% for i = 1:3
%     text(i, final_ipr(i) + 0.02, sprintf('%.4f', final_ipr(i)), ...
%         'HorizontalAlignment', 'center', 'FontName', fontName, 'FontSize', fontSize-2);
% end
% 
% set(gca, 'FontName', fontName, 'FontSize', fontSize);
% set(gca, 'XTickLabel', {'初态1', '初态2', '初态3'});
% xlabel('初态类型', 'FontName', fontName, 'FontSize', fontSize);
% ylabel('最终IPR值', 'FontName', fontName, 'FontSize', fontSize);
% title(sprintf('最终时刻(t=%.1f)的IPR对比', t_max), 'FontName', fontName, 'FontSize', fontSize+2);
% grid on; box on;
% ylim([0, 1]);
% 
% % % 添加理论极限线
% % hold on;
% % plot([0.5, 3.5], [1/total_sites, 1/total_sites], 'k--', 'LineWidth', 1.5, 'DisplayName', sprintf('均匀极限 1/N=%.4f', 1/total_sites));
% % plot([0.5, 3.5], [1, 1], 'r--', 'LineWidth', 1.5, 'DisplayName', '完全局域极限');
% % hold off;
% % 
% % % 添加图例
% % legend({'均匀极限', '完全局域极限'}, 'Location', 'best', 'FontName', fontName, 'FontSize', fontSize-2);
% 
% sgtitle(sprintf('IPR演化对比 (N=%d, 总格点数=%d)', N, total_sites), 'FontSize', 16, 'FontWeight', 'bold');
% print('-dpng','-r600',sprintf('Modelc_IPR_comparison_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));
% 
% fprintf('完成！\n');
toc


% buildHamiltonian_c.m - 构建模型c的哈密顿量
function H = buildHamiltonian_c(N, gamma, J1, J2)
% 输入:
%   N: 每部分的原胞数
%   gamma: 增益损耗强度
%   J1, J2: 跃迁强度
% 输出:
%   H: (4N+1)×(4N+1)的哈密顿量矩阵（稀疏矩阵）

total_sites = 4 * N + 1;  % 4N+1个格点
H = sparse(total_sites, total_sites);  % 预分配稀疏矩阵

% 格点索引:
% 1: a_0
% 2: a_1, 3: b_1
% 4: a_2, 5: b_2
% ...
% 2n: a_n, 2n+1: b_n (对于n=1到2N)

% 1. 第0个格点 (a_0) - 索引1
H(1, 1) = +1i * gamma;  % a_0增益 +iγ

% 2. 第一部分: 前2N个格点 (n=1..N) - 第2到第(2N+1)个格点
for n = 1:N
    a_index = 2*n;      % a_n的索引
    b_index = 2*n+1;    % b_n的索引
    
    % 对角项: 增益和损耗
    H(a_index, a_index) = -1i * gamma;  % a_n损耗 -iγ
    H(b_index, b_index) = +1i * gamma;  % b_n增益 +iγ
    
    % 原胞内跃迁 J1 (a_n <-> b_n)
    H(a_index, b_index) = J1;
    H(b_index, a_index) = J1;
    
    % 原胞间跃迁 J2 (b_n <-> a_{n+1})
    if n < N
        a_next_index = 2*(n+1);  % a_{n+1}
        H(b_index, a_next_index) = J2;
        H(a_next_index, b_index) = J2;
    end
end

% 3. 第二部分: 后2N个格点 (n=N+1..2N) - 第(2N+2)到第(4N+1)个格点
for n = N+1:2*N
    a_index = 2*n;      % a_n的索引
    b_index = 2*n+1;    % b_n的索引
    
    % 对角项: 增益和损耗
    H(a_index, a_index) = +1i * gamma;  % a_n增益 +iγ
    H(b_index, b_index) = -1i * gamma;  % b_n损耗 -iγ
    
    % 原胞内跃迁 J2 (a_n <-> b_n)
    H(a_index, b_index) = J2;
    H(b_index, a_index) = J2;
    
    % 原胞间跃迁 J1 (b_n <-> a_{n+1})
    if n < 2*N
        a_next_index = 2*(n+1);  % a_{n+1}
        H(b_index, a_next_index) = J1;
        H(a_next_index, b_index) = J1;
    end
end

% 4. 特殊连接: a_0与a_1之间的跃迁 J2
H(1, 2) = J2;  % a_0(索引1) <-> a_1(索引2)
H(2, 1) = J2;

% 5. 特殊连接: b_N与a_{N+1}之间的跃迁强度为1
b_N_index = 2*N+1;      % b_N的索引
a_Nplus1_index = 2*N+2;  % a_{N+1}的索引
H(b_N_index, a_Nplus1_index) = 1;
H(a_Nplus1_index, b_N_index) = 1;

end

% solveSpectrum_c.m - 求解本征值能谱（支持大N的稀疏矩阵方法）
function [eigenvalues, eigenvectors] = solveSpectrum_c(H, N)
% 输入:
%   H: 哈密顿量矩阵（稀疏或稠密）
%   N: 原胞数（用于判断是否使用稀疏方法）
% 输出:
%   eigenvalues: 本征值向量
%   eigenvectors: 本征向量矩阵

total_sites = 4 * N + 1;  % 总格点数

if N > 100
    % 对于大N，使用稀疏矩阵特征值求解方法
    fprintf('  使用稀疏方法 (N=%d > 100)...\n', N);
    if total_sites > 2000
        % 对于非常大的矩阵，只求解部分本征值
        num_eig = min(500, total_sites);  % 求解500个本征值
        [eigenvectors, eigenvalues] = eigs(H, num_eig);
    else
        % 对于较大的矩阵，求解全部本征值
        [eigenvectors, eigenvalues] = eig(full(H));
    end
else
    % 对于小N，直接求解全部本征值
    fprintf('  使用直接方法 (N=%d <= 100)...\n', N);
    if issparse(H)
        H = full(H);
    end
    [eigenvectors, eigenvalues] = eig(H);
end

eigenvalues = diag(eigenvalues);
% 按实部排序
[~, idx] = sort(real(eigenvalues));
eigenvalues = eigenvalues(idx);
eigenvectors = eigenvectors(:, idx);

end

% timeEvolve_c.m - 时间演化（通用版本）
function [t, psi_t, prob_density, ipr_evolution] = timeEvolve_c(H, psi0, t_max, dt)
% 输入:
%   H: 哈密顿量矩阵（稀疏或稠密）
%   psi0: 初始态向量
%   t_max: 最大演化时间
%   dt: 时间步长
% 输出:
%   t: 时间向量
%   psi_t: 演化后的态矩阵 (每列对应一个时刻)
%   prob_density: 概率密度矩阵
%   ipr_evolution: IPR随时间演化

if issparse(H)
    H = full(H);  % 矩阵指数需要稠密矩阵
end

% 时间向量
t = 0:dt:t_max;
num_steps = length(t);

% 预分配
psi_t = zeros(length(psi0), num_steps);
prob_density = zeros(length(psi0), num_steps);
ipr_evolution = zeros(1, num_steps);

% 初始态
psi = psi0;
psi_t(:, 1) = psi;
prob_density(:, 1) = abs(psi).^2;

% 计算初始IPR
ipr_evolution(1) = sum(prob_density(:, 1).^2) / (sum(prob_density(:, 1))^2);

% 时间演化 (使用矩阵指数方法)
for i = 2:num_steps
    % 计算 e^{-iHΔt}
    U = expm(-1i * H * dt);
    psi = U * psi;
    
    psi_t(:, i) = psi;
    prob_density(:, i) = abs(psi).^2;
    
    % 计算当前时间步的IPR
    ipr_evolution(i) = sum(prob_density(:, i).^2) / (sum(prob_density(:, i))^2);
end

end
