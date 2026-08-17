% main.m - 主脚本，演示模型a的完整分析流程

%% 参数设置
clear; clc; close all;tic
N = 10;          % 原胞数
gamma = 0;    % 增益损耗强度
J1 = 0.5;       % j-delta
J2 = 1.5;       % j+delta
t_max = 1;     % 演化时间
dt = 0.01;       % 时间步长

%% 1. 构建哈密顿量
fprintf('构建哈密顿量 (N=%d)...\n, time(t=%d)...\n', N, t_max);
H = buildHamiltonian(N, gamma, J1, J2);
total_sites = 4*N + 1;  % 更新：4N+1个格点
if N <= 20
    fprintf('哈密顿量矩阵 (前10×10):\n');
    disp(full(H(1:min(10,total_sites), 1:min(10,total_sites))));
end

%% 2. 求解本征能谱
fprintf('求解本征值...\n');
[real_eigenvalues, eigenvalues, eigenvectors] = solveSpectrum(H, N);

% 绘制改进的能谱图
figure('Name', '本征值能谱', 'Position', [100 100 800 600]);
plot(real_eigenvalues, 'b.', 'MarkerSize', 15);
xlabel('本征态索引'); ylabel('能量实部');
title(sprintf('本征值实部 (N=%d, 总格点数=%d)', N, total_sites));
grid on; box on;axis tight;
print('-dpng','-r600',sprintf('Modela_spectrum_N%d_gamma%.2f_J1%.2f_J2%.2f.png', N, gamma, J1, J2));

% 绘制复平面上的本征值分布
figure('Name', '复平面本征值分布', 'Position', [100 100 600 600]);
plot(real(eigenvalues), imag(eigenvalues), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
xlabel('实部'); ylabel('虚部');
title(sprintf('复平面上的本征值分布 (N=%d, 总格点数=%d)', N, total_sites));
grid on; box on;
axis equal;
print('-dpng','-r600',sprintf('Modela_real_imag_N%d_gamma%.2f_J1%.2f_J2%.2f.png', N, gamma, J1, J2));

%% 3. 时间演化
fprintf('执行时间演化...\n');
psi0 = zeros(total_sites, 1);
b_N_index = 2*N+1;  % b_N的索引（注意：现在有第0个格点，索引需要调整）
psi0(b_N_index) = 1;  % 在b_N上设置初态

% 模拟时间演化
[t, psi_t, prob_b_N, ipr_evolution] = timeEvolve(H, psi0, t_max, dt, b_N_index);

% 绘制概率随时间演化
figure('Name', 'b_N格点概率演化', 'Position', [100 100 800 400]);
plot(t, prob_b_N, 'LineWidth', 2);
xlabel('时间 t'); ylabel('概率 |ψ(b_N)|²');
title(['b_N格点上的概率演化 (N=' num2str(N) ', 总格点数=' num2str(total_sites) ')']);
grid on; box on;
print('-dpng','-r600',sprintf('Modela_prob_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));

% 绘制IPR随时间演化
figure('Name', 'IPR随时间演化', 'Position', [100 100 800 400]);
plot(t, ipr_evolution, 'LineWidth', 2, 'Color', [0.8, 0.2, 0.2]);
xlabel('时间 t'); ylabel('IPR');
title(['逆参与率随时间演化 (N=' num2str(N) ', 总格点数=' num2str(total_sites) ')']);
grid on; box on;
print('-dpng','-r600',sprintf('Modela_IPR_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));

% % 输出IPR统计信息
% fprintf('\nIPR统计信息:\n');
% fprintf('初始IPR: %.6f\n', ipr_evolution(1));
% fprintf('最终IPR: %.6f\n', ipr_evolution(end));
% fprintf('平均IPR: %.6f\n', mean(ipr_evolution));
% fprintf('IPR变化范围: [%.6f, %.6f]\n', min(ipr_evolution), max(ipr_evolution));

%% 4. 可视化演化过程
% 计算概率密度矩阵
prob_density = abs(psi_t).^2;

% 热图：展示所有格点的概率密度随时间演化
figure('Name', '概率密度演化热图', 'Position', [100 100 900 400]);
imagesc(1:total_sites, t, prob_density');
colorbar;
ylabel('时间');
xlabel('格点位置');
title(sprintf('所有格点的概率密度随时间演化 (N=%d, 总格点数=%d)', N, total_sites));
axis xy;
grid on; box on;
print('-dpng','-r600',sprintf('Modela_density_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));

% % 标记特殊格点
% hold on;
% plot(b_N_index, 0, 'w.', 'MarkerSize', 15);  % 第0个格点和b_N
% % text(1, -0.5, '格点0', 'Color', 'w', 'HorizontalAlignment', 'center');
% text(b_N_index, -0.5, 'b_N', 'Color', 'w', 'HorizontalAlignment', 'center');
% hold off;

% 最终时刻的概率分布
figure('Name', '最终概率分布', 'Position', [100 100 800 400]);
bar(1:total_sites, prob_density(:, end), 'FaceColor', [0.5, 0.8, 1.0]);
xlabel('格点位置'); ylabel('概率');
title(sprintf('t=%.1f时的概率分布 (总格点数=%d)', t_max, total_sites));
grid on; box on;
xlim([0, total_sites + 1]);
hold on;
% plot(1, prob_density(1, end), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(b_N_index, prob_density(b_N_index, end), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
% text(1, prob_density(1, end), '  格点0', 'Color', 'g');
text(b_N_index, prob_density(b_N_index, end), '  b_N', 'Color', 'r');
hold off;
print('-dpng','-r600',sprintf('Modela_final_N%d_gamma%.2f_J1%.2f_J2%.2f_t%d.png', N, gamma, J1, J2, t_max));

fprintf('完成！\n');
%% 5. 分析γ=0时的零能本征态
fprintf('\n分析γ=0时的零能本征态...\n');

% 设置γ=0
gamma_zero = 0;

% 构建γ=0时的哈密顿量
H_zero = buildHamiltonian(N, gamma_zero, J1, J2);

% 求解本征值和本征态
[real_eigenvalues_zero, eigenvalues_zero, eigenvectors_zero] = solveSpectrum(H_zero, N);

% 找到实部接近0的本征值索引
zero_threshold = 1e-10;  % 零能阈值
zero_indices = find(abs(real(eigenvalues_zero)) < zero_threshold);

if ~isempty(zero_indices)
    fprintf('找到 %d 个实部接近0的本征态:\n', length(zero_indices));
    
    % 创建零能本征态分布图
    figure('Name', 'γ=0时的零能本征态分布', 'Position', [100 100 1200 800]);
    
    num_zero_states = length(zero_indices);
    rows = ceil(sqrt(num_zero_states));
    cols = ceil(num_zero_states / rows);
    
    for i = 1:num_zero_states
        idx = zero_indices(i);
        eigenstate = eigenvectors_zero(:, idx);
        
        % 计算概率分布
        prob_dist = abs(eigenstate).^2;
        
        % 归一化概率分布
        prob_dist = prob_dist / sum(prob_dist);
        
        % 绘制子图
        subplot(rows, cols, i);
        bar(1:total_sites, prob_dist, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'none');
        
        % 标记特殊格点
        hold on;
        b_N_index = 2*N+1;
        a_Nplus1_index = 2*N+2;
        
        % 标记第0个格点
        plot(1, prob_dist(1), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
        text(1, prob_dist(1)*1.1, 'a_0', 'Color', 'm', ...
            'HorizontalAlignment', 'center', 'FontSize', 10);
        
        % 标记b_N
        if b_N_index <= total_sites
            plot(b_N_index, prob_dist(b_N_index), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
            text(b_N_index, prob_dist(b_N_index)*1.1, 'b_N', 'Color', 'r', ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
        end
        
        % 标记a_{N+1}
        if a_Nplus1_index <= total_sites
            plot(a_Nplus1_index, prob_dist(a_Nplus1_index), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
            text(a_Nplus1_index, prob_dist(a_Nplus1_index)*1.1, 'a_{N+1}', 'Color', 'g', ...
                'HorizontalAlignment', 'center', 'FontSize', 10);
        end
        hold off;
        
        % 计算IPR
        ipr_value = sum(prob_dist.^2);
        
        % 设置子图属性
        xlabel('格点位置', 'FontSize', 10);
        ylabel('概率', 'FontSize', 10);
        title(sprintf('零能态 %d\nIPR=%.4f', i, ipr_value), 'FontSize', 11);
        grid on; box on;
        xlim([0, total_sites+1]);
        ylim([0, max(prob_dist)*1.2]);
    end
    
    sgtitle(sprintf('γ=0时的零能本征态分布 (N=%d, J1=%.2f, J2=%.2f, 总格点数=%d)', N, J1, J2, total_sites), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % 保存图像
    print('-dpng', '-r600', sprintf('Modela_gamma0_zero_states_N%d_J1%.2f_J2%.2f.png', N, J1, J2));
    
    % 创建零能本征态热图对比
    figure('Name', '零能本征态热图对比', 'Position', [100 100 1400 600]);
    
    % 将所有零能态的概率分布组合成矩阵
    zero_states_matrix = zeros(total_sites, num_zero_states);
    for i = 1:num_zero_states
        idx = zero_indices(i);
        eigenstate = eigenvectors_zero(:, idx);
        prob_dist = abs(eigenstate).^2;
        prob_dist = prob_dist / sum(prob_dist);  % 归一化
        zero_states_matrix(:, i) = prob_dist;
    end
    
    % 绘制热图
    imagesc(1:num_zero_states, 1:total_sites, zero_states_matrix);
    colorbar;
    xlabel('零能态索引', 'FontSize', 12);
    ylabel('格点位置', 'FontSize', 12);
    title(sprintf('γ=0时零能本征态的概率分布热图 (N=%d, 总格点数=%d)', N, total_sites), 'FontSize', 14);
    axis xy;
    grid on; box on;
    
    % 标记特殊格点
    hold on;
    plot([0.5, num_zero_states+0.5], [1, 1], 'm-', 'LineWidth', 2);  % a_0
    plot([0.5, num_zero_states+0.5], [b_N_index, b_N_index], 'r-', 'LineWidth', 2);  % b_N
    plot([0.5, num_zero_states+0.5], [a_Nplus1_index, a_Nplus1_index], 'g-', 'LineWidth', 2);  % a_{N+1}
    text(num_zero_states+0.2, 1, 'a_0', 'Color', 'm', 'FontSize', 12, 'FontWeight', 'bold');
    text(num_zero_states+0.2, b_N_index, 'b_N', 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');
    text(num_zero_states+0.2, a_Nplus1_index, 'a_{N+1}', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');
    hold off;
    
    % 保存热图
    print('-dpng', '-r600', sprintf('Modela_gamma0_zero_states_heatmap_N%d_J1%.2f_J2%.2f.png', N, J1, J2));
    
    % 输出零能态的详细信息
    fprintf('\n零能本征态详细信息:\n');
    fprintf('%-10s %-15s %-15s %-10s\n', '索引', '能量实部', '能量虚部', 'IPR');
    fprintf('%-10s %-15s %-15s %-10s\n', '------', '----------', '----------', '------');
    
    for i = 1:num_zero_states
        idx = zero_indices(i);
        eigenstate = eigenvectors_zero(:, idx);
        prob_dist = abs(eigenstate).^2;
        prob_dist = prob_dist / sum(prob_dist);
        ipr_value = sum(prob_dist.^2);
        
        fprintf('%-10d %-15.6e %-15.6e %-10.6f\n', ...
            idx, real(eigenvalues_zero(idx)), imag(eigenvalues_zero(idx)), ipr_value);
    end
    
else
    fprintf('未找到实部接近0的本征态。\n');
    fprintf('最接近0的本征值实部: %.6e\n', min(abs(real(eigenvalues_zero))));
end

fprintf('完成！\n');
toc


% buildHamiltonian.m - 构建模型a的哈密顿量（包含第0个格点）
function H = buildHamiltonian(N, gamma, J1, J2)
% 输入:
%   N: 每部分的原胞数
%   gamma: 增益损耗强度
%   J1, J2: 跃迁强度
% 输出:
%   H: (4N+1)×(4N+1)的哈密顿量矩阵（稀疏矩阵）

total_sites = 4 * N + 1;  % 更新：4N+1个格点
H = sparse(total_sites, total_sites);  % 预分配稀疏矩阵

% 第0个格点（索引1）
H(1, 1) = +1i * gamma;  % 第0个格点增益 +iγ

% 第一部分: 前2N个格点 (n=1..N) - 注意：索引从2开始
for n = 1:N
    a_index = 2*n;      % a_n的索引（注意：+1因为第0个格点）
    b_index = 2*n+1;  % b_n的索引
    
    % 对角项: 增益和损耗
    if n < N
        H(a_index, a_index) = -1i * gamma;  % a_n损耗
        H(b_index, b_index) = +1i * gamma;  % b_n增益
    else  % n == N, 特殊处理b_N
        H(a_index, a_index) = -1i * gamma;  % a_N损耗
        H(b_index, b_index) = 0;            % b_N增益损耗为0
    end
    
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

% 第0个格点与a_1之间的跃迁 J2
H(1, 2) = J2;  % 第0个格点(索引1) <-> a_1(索引2)
H(2, 1) = J2;

% 第二部分: 后2N个格点 (n=N+1..2N)
for n = N+1:2*N
    a_index = 2*n;      % a_n的索引
    b_index = 2*n + 1;  % b_n的索引
    
    % 对角项: 增益和损耗
    H(a_index, a_index) = -1i * gamma;  % a_n增益
    H(b_index, b_index) = +1i * gamma;  % b_n损耗
    
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

% 处理两部分连接: b_N (索引2N+1) 与 a_{N+1} (索引2N+2)
H(2*N+1, 2*N+2) = J1;
H(2*N+2, 2*N+1) = J1;

end

% solveSpectrum.m - 求解本征值能谱
function [real_eigenvalues, eigenvalues, eigenvectors] = solveSpectrum(H, N)
% 输入:
%   H: 哈密顿量矩阵（稀疏或稠密）
%   N: 原胞数（用于判断是否使用稀疏方法）
% 输出:
%   eigenvalues: 本征值向量
%   eigenvectors: 本征向量矩阵

total_sites = 4 * N + 1;  % 更新：4N+1个格点

if N > 100
    % 对于大N，转换为稠密矩阵求解（稀疏特征值问题需要特殊处理）
    fprintf('  使用稀疏方法 (N=%d > 100)...\n', N);
    H_dense = full(H);  % 转换为稠密矩阵
    [eigenvectors, eigenvalues] = eig(H_dense);  % 注意: eig返回[D,V] = eig(A)
else
    % 对于小N，直接求解
    fprintf('  使用直接方法 (N=%d <= 100)...\n', N);
    if issparse(H)
        H = full(H);
    end
    [eigenvectors, eigenvalues] = eig(H);
end

eigenvalues = diag(eigenvalues);
% 按实部排序
[real_eigenvalues,~]=sort(real(eigenvalues));

end

% timeEvolve.m - 时间演化
function [t, psi_t, prob_snapshot, ipr_evolution] = timeEvolve(H, psi0, t_max, dt, target_index)
% 输入:
%   H: 哈密顿量矩阵（稀疏或稠密）
%   psi0: 初始态向量
%   t_max: 最大演化时间
%   dt: 时间步长
%   target_index: 要监测的格点索引
% 输出:
%   t: 时间向量
%   psi_t: 演化后的态矩阵 (每列对应一个时刻)
%   prob_snapshot: 目标格点的概率随时间演化
%   ipr_evolution: IPR随时间演化

if issparse(H)
    H = full(H);  % 矩阵指数需要稠密矩阵
end

% 时间向量
t = 0:dt:t_max;
num_steps = length(t);

% 预分配
psi_t = zeros(length(psi0), num_steps);
prob_snapshot = zeros(1, num_steps);
ipr_evolution = zeros(1, num_steps);

% 初始态
psi = psi0;
psi_t(:, 1) = psi;
prob_snapshot(1) = abs(psi(target_index))^2;

% 计算初始IPR
prob_density = abs(psi).^2;
ipr_evolution(1) = sum(prob_density.^2) / (sum(prob_density)^2);

% 时间演化 (使用矩阵指数方法)
for i = 2:num_steps
    % 计算 e^{-iHΔt}
    U = expm(-1i * H * dt);
    psi = U * psi;
    
    psi_t(:, i) = psi;
    prob_snapshot(i) = abs(psi(target_index))^2;
    
    % 计算当前时间步的IPR
    prob_density = abs(psi).^2;
    ipr_evolution(i) = sum(prob_density.^2) / (sum(prob_density)^2);
end

end