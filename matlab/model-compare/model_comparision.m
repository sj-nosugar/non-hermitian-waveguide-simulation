% 比较模型a、b和c的末态分布及演化热图
% 模型c包含从第N个、第N+1个格点开始及两者各0.5的初始状态
clear; clc; close all; tic    
    % -------------------------- 修正：确保文件夹可识别 --------------------------
    fig_folder = 'result-fig';  % 文件夹名称
    % 检查文件夹是否存在，不存在则创建
    if ~exist(fig_folder, 'dir')
        mkdir(fig_folder);
        fprintf('已创建图片保存文件夹：%s\n', fullfile(pwd, fig_folder));
    else
        fprintf('图片保存文件夹已存在：%s\n', fullfile(pwd, fig_folder));
    end

    % 统一参数设置
    N = 4;          % 格点参数
    t = 1;           % 跃迁强度
    delta = 0.5;     % 跃迁扰动
    gamma = 0;     % 增益耗散强度
    total_time = 1; % 总演化时间
    num_steps = 100; % 时间步数
    
    % 运行模型a
    fprintf('正在运行模型a...\n');
    [~, ~, ~] = construct_hamiltonian_a(N, t, delta, gamma);
    [~, prob_a] = simulate_evolution_a(N, t, delta, gamma, total_time, num_steps);
    
    % 运行模型b
    fprintf('正在运行模型b...\n');
    [~, ~, ~] = construct_hamiltonian_b(N, t, delta, gamma);
    [~, prob_b] = simulate_evolution_b(N, t, delta, gamma, total_time, num_steps);
    
    % 运行模型c（三种初始状态）
    fprintf('正在运行模型c...\n');
    [~, ~, ~] = construct_hamiltonian_c(N, t, delta, gamma);
    [~, ~, prob_c_N, prob_c_N1, prob_c_N_N1] = simulate_evolution_c(N, t, delta, gamma, total_time, num_steps);
    
    % 生成对比可视化
    fprintf('正在生成对比可视化结果...\n');
    plot_final_distributions(prob_a, prob_b, prob_c_N, prob_c_N1, prob_c_N_N1, N, fig_folder);
    plot_evolution_heatmaps(prob_a, prob_b, prob_c_N, prob_c_N1, prob_c_N_N1, N, total_time, num_steps, fig_folder);
    
    fprintf('对比分析完成！所有图片已保存至：%s\n', fullfile(pwd, fig_folder));
toc

%%-------- 末态分布绘图函数 --------
function plot_final_distributions(prob_a, prob_b, prob_c_N, prob_c_N1, prob_c_N_N1, N, fig_folder)
    % 绘制所有模型的末态分布对比图
    
    % 获取各模型的总格点数
    total_sites_a_b = 2*N + 1;  % 模型a和b的总格点数
    total_sites_c = 2*N;        % 模型c的总格点数
    sites_a_b = 1:total_sites_a_b;
    sites_c = 1:total_sites_c;
    
    % 1. 综合对比图
    figure('Name', '各模型末态分布对比', 'Position', [100 100 1400 1000]);
    set(gcf,'Color','w','PaperPositionMode','auto');
    fontName = 'Helvetica';
    fontSize = 16;
    
    % 模型a末态分布
    subplot(3, 2, 1);
    bar(sites_a_b, prob_a(end, :), 'FaceColor', [0.5, 0.8, 1.0]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title('模型a末态概率分布','FontSize',fontSize,'FontName',fontName);
    grid on;box on;xlim([0, total_sites_a_b + 1]);ylim([0, 0.7]);
    
    % 模型b末态分布
    subplot(3, 2, 2);
    bar(sites_a_b, prob_b(end, :), 'FaceColor', [0.8, 0.5, 0.8]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title('模型b末态概率分布','FontSize',fontSize,'FontName',fontName);
    grid on;box on;xlim([0, total_sites_a_b + 1]);ylim([0, 0.7]);
    
    % 模型c（从N开始）末态分布
    subplot(3, 2, 3);
    bar(sites_c, prob_c_N(end, :), 'FaceColor', [0.5, 0.8, 0.5]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title('模型c（从第N个格点开始）','FontSize',fontSize,'FontName',fontName);
    grid on;box on;xlim([0, total_sites_c + 1]);ylim([0, 0.7]);
    
    % 模型c（从N+1开始）末态分布
    subplot(3, 2, 4);
    bar(sites_c, prob_c_N1(end, :), 'FaceColor', [1.0, 0.6, 0.5]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title('模型c（从第N+1个格点开始）','FontSize',fontSize,'FontName',fontName);
    grid on;box on;xlim([0, total_sites_c + 1]);ylim([0, 0.7]);
    
    % 模型c（N和N+1各0.5）末态分布
    subplot(3, 2, 5);
    bar(sites_c, prob_c_N_N1(end, :), 'FaceColor', [0.9, 0.9, 0.2]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title('模型c（N和N+1各0.5）','FontSize',fontSize,'FontName',fontName);
    grid on;box on;xlim([0, total_sites_c + 1]);ylim([0, 0.7]);
    
    % 导出高清TIFF到result-fig文件夹
    save_path = fullfile(fig_folder, 'Model-final-comparison.tiff');
    print(gcf, '-dtiff', '-r600', save_path);
    fprintf('末态分布综合对比图已保存：%s\n', save_path);
end

%%-------- 演化热图绘图函数 --------
function plot_evolution_heatmaps(prob_a, prob_b, prob_c_N, prob_c_N1, prob_c_N_N1, N, total_time, num_steps, fig_folder)
    % 绘制所有模型的演化热图对比
    
    % 时间点与格点位置
    time_points = linspace(0, total_time, num_steps + 1);
    sites_a_b = 1:1:(2*N + 1);
    sites_c = 1:1:(2*N);
    fontName = 'Helvetica';
    fontSize = 16;
    
    % 1. 综合热图对比
    figure('Name', '各模型演化热图对比', 'Position', [100 100 1400 1000]);
    set(gcf,'Color','w','PaperPositionMode','auto');
    
    % 模型a热图
    subplot(3, 2, 1);
    imagesc(sites_a_b, time_points, prob_a);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;ylabel('时间','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title('模型a概率密度演化','FontSize',fontSize,'FontName',fontName);
    axis xy;grid on;box on;clim([0, 1]);
    
    % 模型b热图
    subplot(3, 2, 2);
    imagesc(sites_a_b, time_points, prob_b);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;ylabel('时间','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title('模型b概率密度演化','FontSize',fontSize,'FontName',fontName);
    axis xy;grid on;box on;clim([0, 1]);
    
    % 模型c（从N开始）热图
    subplot(3, 2, 3);
    imagesc(sites_c, time_points, prob_c_N);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;ylabel('时间','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title('模型c（从第N个格点开始）','FontSize',fontSize,'FontName',fontName);
    axis xy;grid on;box on;clim([0, 1]);
    
    % 模型c（从N+1开始）热图
    subplot(3, 2, 4);
    imagesc(sites_c, time_points, prob_c_N1);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;ylabel('时间','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title('模型c（从第N+1个格点开始）','FontSize',fontSize,'FontName',fontName);
    axis xy;grid on;box on;clim([0, 1]);
    
    % 模型c（N和N+1各0.5）热图
    subplot(3, 2, 5);
    imagesc(sites_c, time_points, prob_c_N_N1);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;ylabel('时间','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title('模型c（N和N+1各0.5）','FontSize',fontSize,'FontName',fontName);
    axis xy;grid on;box on;clim([0, 1]);
    
    % 导出热图
    save_path = fullfile(fig_folder, 'Models_Evolution-comparison.tiff');
    print(gcf, '-dtiff', '-r600', save_path);
    fprintf('演化热图综合对比图已保存：%s\n', save_path);
end

% =================== 模型a函数 ===================
function [H, eigenvalues, eigenvectors] = construct_hamiltonian_a(N, t, delta, gamma)
    % 构造模型a的哈密顿量并求解本征值和本征向量
    total_sites = 2*N + 1;
    H = sparse(total_sites, total_sites);
    
    for i = 1:total_sites
        if i == N + 1
            H(i, i) = 0;
        else
            if mod(N, 2) == 1  
                if i <= N  
                    if mod(i, 2) == 1  
                        H(i, i) = -1i * gamma;
                    else  
                        H(i, i) = 1i * gamma;
                    end
                else  
                    if mod(i, 2) == 1  
                        H(i, i) = 1i * gamma;
                    else  
                        H(i, i) = -1i * gamma;
                    end
                end
            else  
                if i <= N  
                    if mod(i, 2) == 1  
                        H(i, i) = 1i * gamma;
                    else  
                        H(i, i) = -1i * gamma;
                    end
                else  
                    if mod(i, 2) == 1  
                        H(i, i) = 1i * gamma;
                    else  
                        H(i, i) = -1i * gamma;
                    end
                end
            end
        end
    end
    
    if mod(N, 2) == 1  
        H(1, 2) = t - delta;
        H(2, 1) = conj(t - delta);
        num_cells = (N - 1) / 2;
        for n = 1:num_cells
            site1 = 2*n;
            site2 = 2*n + 1;
            if site2 <= N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            site1 = 2*n + 1;
            site2 = 2*n + 2;
            if site2 <= N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
        H(N, N+1) = t - delta;
        H(N+1, N) = conj(t - delta);
        for n = 1:num_cells
            site1 = N + 1 + 2*n - 2;
            site2 = N + 1 + 2*n - 1;
            if site2 <= total_sites
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
            site1 = N + 1 + 2*n - 1;
            site2 = N + 1 + 2*n;
            if site2 <= total_sites
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
        end
        H(total_sites - 1, total_sites) = t - delta;
        H(total_sites, total_sites - 1) = conj(t - delta);
    else  
        num_cells = N / 2;
        for n = 1:num_cells
            site1 = 2*n - 1;
            site2 = 2*n;
            if site2 <= N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            site1 = 2*n;
            site2 = 2*n + 1;
            if site2 <= N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
        H(N, N+1) = t - delta;
        H(N+1, N) = conj(t - delta);
        for n = 1:num_cells
            site1 = N + 1 + 2*n - 2;
            site2 = N + 1 + 2*n - 1;
            if site1 >= N + 1 && site2 <= total_sites
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
            site1 = N + 1 + 2*n - 1;
            site2 = N + 1 + 2*n;
            if site2 <= total_sites
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
        end
    end
    
    if total_sites <= 1000  
        H_full = full(H);
        [eigenvectors, eigenvalues_matrix] = eig(H_full);
        eigenvalues = diag(eigenvalues_matrix);
    else  
        [eigenvectors, eigenvalues_matrix] = eigs(H, total_sites);
        eigenvalues = diag(eigenvalues_matrix);
    end
end

function [evolution, prob_density] = simulate_evolution_a(N, t, delta, gamma, total_time, num_steps)
    % 模拟模型a的时间演化
    [H, ~, ~] = construct_hamiltonian_a(N, t, delta, gamma);
    total_sites = 2*N + 1;
    dt = total_time / num_steps;
    
    initial_state = zeros(total_sites, 1);
    initial_state(N + 1) = 1;
    
    evolution = zeros(num_steps + 1, total_sites);
    evolution(1, :) = initial_state';
    prob_density = zeros(num_steps + 1, total_sites);
    prob_density(1, :) = abs(initial_state').^2;
    
    current_state = initial_state;
    for step = 1:num_steps
        if total_sites <= 1000
            U = expm(-1i * full(H) * dt);
        else
            U = speye(total_sites) - 1i * H * dt - 0.5 * (H^2) * (dt^2);
        end
        current_state = U * current_state;
        current_state = current_state / norm(current_state);
        evolution(step + 1, :) = current_state';
        prob_density(step + 1, :) = abs(current_state').^2;
    end
end

% =================== 模型b函数 ===================
function [H, eigenvalues, eigenvectors] = construct_hamiltonian_b(N, t, delta, gamma)
    % 构造模型b的哈密顿量并求解本征值和本征向量
    total_sites = 2*N + 1;
    H = sparse(total_sites, total_sites);
    
    for i = 1:total_sites
        if i == N + 1
            H(i, i) = 0;
        else
            if mod(N, 2) == 1  
                if mod(i, 2) == 1  
                    H(i, i) = -1i * gamma;
                else  
                    H(i, i) = 1i * gamma;
                end
            else  
                if i <= N  
                    if mod(i, 2) == 1  
                        H(i, i) = 1i * gamma;
                    else  
                        H(i, i) = -1i * gamma;
                    end
                else  
                    if mod(i, 2) == 1  
                        H(i, i) = -1i * gamma;
                    else  
                        H(i, i) = 1i * gamma;
                    end
                end
            end
        end
    end
    
    if mod(N, 2) == 1  
        H(1, 2) = t - delta;
        H(2, 1) = conj(t - delta);
        num_cells = (N - 1) / 2;
        for n = 1:num_cells
            site1 = 2*n;
            site2 = 2*n + 1;
            if site2 <= N
                H(site1, site2) = t - delta;  % 与模型a不同
                H(site2, site1) = conj(t - delta);
            end
            site1 = 2*n + 1;
            site2 = 2*n + 2;
            if site2 <= N
                H(site1, site2) = t + delta;  % 与模型a不同
                H(site2, site1) = conj(t + delta);
            end
        end
        H(N, N+1) = t + delta;  % 与模型a不同
        H(N+1, N) = conj(t + delta);
        for n = 1:num_cells
            site1 = N + 1 + 2*n - 2;
            site2 = N + 1 + 2*n - 1;
            if site2 <= total_sites
                H(site1, site2) = t + delta;  % 与模型a不同
                H(site2, site1) = conj(t + delta);
            end
            site1 = N + 1 + 2*n - 1;
            site2 = N + 1 + 2*n;
            if site2 <= total_sites
                H(site1, site2) = t - delta;  % 与模型a不同
                H(site2, site1) = conj(t - delta);
            end
        end
        H(total_sites - 1, total_sites) = t + delta;  % 与模型a不同
        H(total_sites, total_sites - 1) = conj(t + delta);
    else  
        num_cells = N / 2;
        for n = 1:num_cells
            site1 = 2*n - 1;
            site2 = 2*n;
            if site2 <= N
                H(site1, site2) = t - delta;  % 与模型a不同
                H(site2, site1) = conj(t - delta);
            end
            site1 = 2*n;
            site2 = 2*n + 1;
            if site2 <= N
                H(site1, site2) = t + delta;  % 与模型a不同
                H(site2, site1) = conj(t + delta);
            end
        end
        H(N, N+1) = t + delta;  % 与模型a不同
        H(N+1, N) = conj(t + delta);
        for n = 1:num_cells
            site1 = N + 1 + 2*n - 2;
            site2 = N + 1 + 2*n - 1;
            if site1 >= N + 1 && site2 <= total_sites
                H(site1, site2) = t + delta;  % 与模型a不同
                H(site2, site1) = conj(t + delta);
            end
            site1 = N + 1 + 2*n - 1;
            site2 = N + 1 + 2*n;
            if site2 <= total_sites
                H(site1, site2) = t - delta;  % 与模型a不同
                H(site2, site1) = conj(t - delta);
            end
        end
    end
    
    if total_sites <= 1000  
        H_full = full(H);
        [eigenvectors, eigenvalues_matrix] = eig(H_full);
        eigenvalues = diag(eigenvalues_matrix);
    else  
        [eigenvectors, eigenvalues_matrix] = eigs(H, total_sites);
        eigenvalues = diag(eigenvalues_matrix);
    end
end

function [evolution, prob_density] = simulate_evolution_b(N, t, delta, gamma, total_time, num_steps)
    % 模拟模型b的时间演化
    [H, ~, ~] = construct_hamiltonian_b(N, t, delta, gamma);
    total_sites = 2*N + 1;
    dt = total_time / num_steps;
    
    initial_state = zeros(total_sites, 1);
    initial_state(N + 1) = 1;
    
    evolution = zeros(num_steps + 1, total_sites);
    evolution(1, :) = initial_state';
    prob_density = zeros(num_steps + 1, total_sites);
    prob_density(1, :) = abs(initial_state').^2;
    
    current_state = initial_state;
    for step = 1:num_steps
        if total_sites <= 1000
            U = expm(-1i * full(H) * dt);
        else
            U = speye(total_sites) - 1i * H * dt - 0.5 * (H^2) * (dt^2);
        end
        current_state = U * current_state;
        current_state = current_state / norm(current_state);
        evolution(step + 1, :) = current_state';
        prob_density(step + 1, :) = abs(current_state').^2;
    end
end

% =================== 模型c函数 ===================
function [H, eigenvalues, eigenvectors] = construct_hamiltonian_c(N, t, delta, gamma)
    % 构造模型c的哈密顿量并求解本征值和本征向量
    total_sites = 2*N;
    H = sparse(total_sites, total_sites);
    
    for i = 1:total_sites
        if i <= N  
            if mod(N, 2) == 1  
                if mod(i, 2) == 1  
                    H(i, i) = 1i * gamma;
                else  
                    H(i, i) = -1i * gamma;
                end
            else  
                if mod(i, 2) == 1  
                    H(i, i) = -1i * gamma;
                else  
                    H(i, i) = 1i * gamma;
                end
            end
        else  
            if mod(i, 2) == 1  
                H(i, i) = 1i * gamma;
            else  
                H(i, i) = -1i * gamma;
            end
        end
    end
    
    H(N, N+1) = t;
    H(N+1, N) = conj(t);
    
    if mod(N, 2) == 1  
        num_cells = (N - 1) / 2;
        for cell = 1:num_cells
            site1 = 2*cell - 1;
            site2 = 2*cell;
            if site2 <= N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            site1 = 2*cell;
            site2 = 2*cell + 1;
            if site2 <= N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
        H(2*N-1, 2*N) = t - delta;
        H(2*N, 2*N-1) = conj(t - delta);
        for cell = 1:num_cells
            site1 = N + 1 + 2*(cell - 1);
            site2 = N + 1 + 2*(cell - 1) + 1;
            if site2 <= 2*N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            site1 = N + 1 + 2*(cell - 1) + 1;
            site2 = N + 1 + 2*(cell - 1) + 2;
            if site2 <= 2*N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
    else  
        H(1, 2) = t - delta;
        H(2, 1) = conj(t - delta);
        H(N-1, N) = t - delta;
        H(N, N-1) = conj(t - delta);
        if N >= 2
            num_cells = (N - 2) / 2;
            for cell = 1:num_cells
                site1 = 2 + 2*(cell - 1);
                site2 = 2 + 2*(cell - 1) + 1;
                if site2 <= N
                    H(site1, site2) = t + delta;
                    H(site2, site1) = conj(t + delta);
                end
                site1 = 2 + 2*(cell - 1) + 1;
                site2 = 2 + 2*(cell - 1) + 2;
                if site2 <= N
                    H(site1, site2) = t - delta;
                    H(site2, site1) = conj(t - delta);
                end
            end
        end
        num_cells = N / 2;
        for cell = 1:num_cells
            site1 = N + 1 + 2*(cell - 1);
            site2 = N + 1 + 2*(cell - 1) + 1;
            if site2 <= 2*N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            site1 = N + 1 + 2*(cell - 1) + 1;
            site2 = N + 1 + 2*(cell - 1) + 2;
            if site2 <= 2*N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
    end
    
    if total_sites <= 1000  
        H_full = full(H);
        [eigenvectors, eigenvalues_matrix] = eig(H_full);
        eigenvalues = diag(eigenvalues_matrix);
    else  
        [eigenvectors, eigenvalues_matrix] = eigs(H, total_sites);
        eigenvalues = diag(eigenvalues_matrix);
    end
end

function [evolution_N, evolution_N1, prob_N, prob_N1, prob_N_N1] = simulate_evolution_c(N, t, delta, gamma, total_time, num_steps)
    % 模拟模型c的时间演化（三种初始状态）
    [H, ~, ~] = construct_hamiltonian_c(N, t, delta, gamma);
    total_sites = 2*N;
    dt = total_time / num_steps;
    
    % 初始状态1：第N个格点
    initial_state_N = zeros(total_sites, 1);
    initial_state_N(N) = 1;
    
    % 初始状态2：第N+1个格点
    initial_state_N1 = zeros(total_sites, 1);
    initial_state_N1(N+1) = 1;
    
    % 新增初始状态3：第N和N+1个格点各0.5概率（已归一化）
    initial_state_N_N1 = zeros(total_sites, 1);
    initial_state_N_N1(N) = 1/sqrt(2);   % 振幅为1/√2，概率为0.5
    initial_state_N_N1(N+1) = 1/sqrt(2); % 振幅为1/√2，概率为0.5
    
    % 初始化存储变量
    evolution_N = zeros(num_steps + 1, total_sites);
    evolution_N(1, :) = initial_state_N';
    prob_N = zeros(num_steps + 1, total_sites);
    prob_N(1, :) = abs(initial_state_N').^2;
    
    evolution_N1 = zeros(num_steps + 1, total_sites);
    evolution_N1(1, :) = initial_state_N1';
    prob_N1 = zeros(num_steps + 1, total_sites);
    prob_N1(1, :) = abs(initial_state_N1').^2;
    
    % 新增：第三种初态的存储变量
    evolution_N_N1 = zeros(num_steps + 1, total_sites);
    evolution_N_N1(1, :) = initial_state_N_N1';
    prob_N_N1 = zeros(num_steps + 1, total_sites);
    prob_N_N1(1, :) = abs(initial_state_N_N1').^2;
    
    % 计算演化算符
    if total_sites <= 1000
        U = expm(-1i * full(H) * dt);
    else
        U = speye(total_sites) - 1i * H * dt - 0.5 * (H^2) * (dt^2);
    end
    
    % 三种初态同时演化
    current_state_N = initial_state_N;
    current_state_N1 = initial_state_N1;
    current_state_N_N1 = initial_state_N_N1;
    
    for step = 1:num_steps
        current_state_N = U * current_state_N;
        current_state_N1 = U * current_state_N1;
        current_state_N_N1 = U * current_state_N_N1;
        
        % 归一化
        current_state_N = current_state_N / norm(current_state_N);
        current_state_N1 = current_state_N1 / norm(current_state_N1);
        current_state_N_N1 = current_state_N_N1 / norm(current_state_N_N1);
        
        % 存储结果
        evolution_N(step + 1, :) = current_state_N';
        prob_N(step + 1, :) = abs(current_state_N').^2;
        
        evolution_N1(step + 1, :) = current_state_N1';
        prob_N1(step + 1, :) = abs(current_state_N1').^2;
        
        % 新增：第三种初态的结果存储
        evolution_N_N1(step + 1, :) = current_state_N_N1';
        prob_N_N1(step + 1, :) = abs(current_state_N_N1').^2;
    end
end

function w = expmv(t, A, v, tol)
    % 使用Krylov子空间方法计算 exp(t*A)*v
    if nargin < 4
        tol = 1e-10;
    end
    
    n = length(v);
    m = min(30, n);
    beta = norm(v);
    v = v / beta;
    
    V = zeros(n, m+1);
    H = zeros(m+1, m);
    V(:, 1) = v;
    
    for j = 1:m
        w = A * V(:, j);
        for i = 1:j
            H(i, j) = V(:, i)' * w;
            w = w - H(i, j) * V(:, i);
        end
        H(j+1, j) = norm(w);
        if H(j+1, j) < tol
            m = j;
            break;
        end
        V(:, j+1) = w / H(j+1, j);
    end
    
    e1 = zeros(m, 1);
    e1(1) = 1;
    expH_e1 = expm(t * H(1:m, 1:m)) * e1;
    w = beta * V(:, 1:m) * expH_e1;
end
