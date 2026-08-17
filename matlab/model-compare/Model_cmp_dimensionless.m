% 比较模型a、b和c的末态分布及演化热图
    % 模型c包含从第N个和第N+1个格点开始的两种初始状态
clear;clc;close all;toc    
    % 统一参数设置
    N = 4;          % 格点参数
    t = 1;           % 跃迁强度
    delta = 0.5;     % 跃迁扰动
    gamma = 0;     % 增益耗散强度
    total_time = 1; % 总演化时间
    num_steps = 100; % 时间步数
    
    % 运行模型a 
    fprintf('正在运行模型a...\n');
    [H_a, eigenvalues_a, ~] = construct_hamiltonian_a(N, t, delta, gamma);
    E_scale_a = max(abs(eigenvalues_a)) - min(abs(eigenvalues_a)); % 计算能量尺度
    if E_scale_a < 1e-10, E_scale_a = 1; end
    t_max_a = total_time / E_scale_a; % 调整真实时间使无量纲时间一致
    [~, prob_a] = simulate_evolution_a(N, t, delta, gamma, t_max_a, num_steps);
    
    % 运行模型b
    fprintf('正在运行模型b...\n');
    [H_b, eigenvalues_b, ~] = construct_hamiltonian_b(N, t, delta, gamma);
    E_scale_b = max(abs(eigenvalues_b)) - min(abs(eigenvalues_b));
    if E_scale_b < 1e-10, E_scale_b = 1; end
    t_max_b = total_time / E_scale_b;
    [~, prob_b] = simulate_evolution_b(N, t, delta, gamma, t_max_b, num_steps);
    
    % 运行模型c（三种初始状态）
    fprintf('正在运行模型c...\n');
    [H_c, eigenvalues_c, ~] = construct_hamiltonian_c(N, t, delta, gamma);
    E_scale_c = max(abs(eigenvalues_c)) - min(abs(eigenvalues_c));
    if E_scale_c < 1e-10, E_scale_c = 1; end
    t_max_c = total_time / E_scale_c;
    [~, ~, prob_c_N, prob_c_N1, prob_c_NN1] = simulate_evolution_c(N, t, delta, gamma, t_max_c, num_steps);
    
    % 生成对比可视化
    fprintf('正在生成对比可视化结果...\n');
    
    % 1. 末态模式分布对比
    plot_final_distributions(prob_a, prob_b, prob_c_N, prob_c_N1, prob_c_NN1, N, ...
        E_scale_a, E_scale_b, E_scale_c, total_time);
    
    % 2. 演化热图对比
    plot_evolution_heatmaps(prob_a, prob_b, prob_c_N, prob_c_N1, prob_c_NN1, N, ...
        t_max_a, t_max_b, t_max_c, num_steps, E_scale_a, E_scale_b, E_scale_c, total_time);
    
    fprintf('对比分析完成！\n');
toc
%%-------- function--------
function plot_final_distributions(prob_a, prob_b, prob_c_N, prob_c_N1, prob_c_NN1, N, ...
    E_scale_a, E_scale_b, E_scale_c, total_time)
    % 绘制所有模型的末态分布对比图
    
    % 获取各模型的总格点数
    total_sites_a_b = 2*N + 1;  % 模型a和b的总格点数
    total_sites_c = 2*N;        % 模型c的总格点数
    sites_a_b = 1:total_sites_a_b;
    sites_c = 1:total_sites_c;
    
    % 创建对比图
    figure('Name', '各模型末态分布对比', 'Position', [100 100 1400 1000]);
    set(gcf,'Color','w','PaperPositionMode','auto');
    fontName = 'Helvetica';
    fontSize = 14;
    
    % 模型a末态分布
    subplot(2, 3, 1);
    bar(sites_a_b, prob_a(end, :), 'FaceColor', [0.5, 0.8, 1.0]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型a末态分布 (E_{scale}=%.3f)', E_scale_a),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_a_b + 1]);
    ylim([0, 0.7]);
    
    % 模型b末态分布
    subplot(2, 3, 2);
    bar(sites_a_b, prob_b(end, :), 'FaceColor', [0.8, 0.5, 0.8]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型b末态分布 (E_{scale}=%.3f)', E_scale_b),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_a_b + 1]);
    ylim([0, 0.7]);
    
    % 模型c（从N开始）末态分布
    subplot(2, 3, 3);
    bar(sites_c, prob_c_N(end, :), 'FaceColor', [0.5, 0.8, 0.5]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（从第N个格点开始）(E_{scale}=%.3f)', E_scale_c),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_c + 1]);
    ylim([0, 0.7]);
    
    % 模型c（从N+1开始）末态分布
    subplot(2, 3, 4);
    bar(sites_c, prob_c_N1(end, :), 'FaceColor', [1.0, 0.6, 0.5]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（从第N+1个格点开始）(E_{scale}=%.3f)', E_scale_c),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_c + 1]);
    ylim([0, 0.7]);
    
    % 模型c（从N和N+1各0.5开始）末态分布
    subplot(2, 3, 5);
    bar(sites_c, prob_c_NN1(end, :), 'FaceColor', [0.9, 0.8, 0.3]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（N和N+1各0.5）(E_{scale}=%.3f)', E_scale_c),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_c + 1]);
    ylim([0, 0.7]);
    
    % 添加总标题
    sgtitle(sprintf('各模型末态分布对比 (无量纲时间Et=%.1f)', total_time), 'FontSize', 16, 'FontWeight', 'bold');
    
    print('-dtiff','-r600','Model-final.tiff');
    
    % 单独图像输出（保持原代码结构）
    figure('Name', 'Model-a', 'Position', [100 100 1200 800]);
    bar(sites_a_b, prob_a(end, :), 'FaceColor', [0.5, 0.8, 1.0]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型a末态概率分布 (Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_a_b + 1]);
    ylim([0, 0.7]);
    print('-dtiff','-r600','Model-a-final.tiff');

    figure('Name', 'Model-b', 'Position', [100 100 1200 800]);
    bar(sites_a_b, prob_b(end, :), 'FaceColor', [0.8, 0.5, 0.8]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型b末态概率分布 (Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_a_b + 1]);
    ylim([0, 0.7]);
    print('-dtiff','-r600','Model-b-final.tiff');
    
    figure('Name', 'Model-c-all', 'Position', [100 100 1800 600]);
    subplot(1,3,1)
    bar(sites_c, prob_c_N(end, :), 'FaceColor', [0.5, 0.8, 0.5]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（从第N个格点开始）(Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_c + 1]);
    ylim([0, 0.7]);
    
    subplot(1, 3, 2);
    bar(sites_c, prob_c_N1(end, :), 'FaceColor', [1.0, 0.6, 0.5]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（从第N+1个格点开始）(Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_c + 1]);
    ylim([0, 0.7]);
    
    subplot(1, 3, 3);
    bar(sites_c, prob_c_NN1(end, :), 'FaceColor', [0.9, 0.8, 0.3]);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    ylabel('概率','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（N和N+1各0.5）(Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    grid on;
    box on;
    xlim([0, total_sites_c + 1]);
    ylim([0, 0.7]);
    print('-dtiff','-r600','Model-c-final.tiff');
end

function plot_evolution_heatmaps(prob_a, prob_b, prob_c_N, prob_c_N1, prob_c_NN1, N, ...
    t_max_a, t_max_b, t_max_c, num_steps, E_scale_a, E_scale_b, E_scale_c, total_time)
    % 绘制所有模型的演化热图对比
    
    % 创建无量纲时间点（所有模型统一为0到total_time）
    time_points = linspace(0, total_time, num_steps + 1);
    
    % 格点位置
    sites_a_b = 1:1:(2*N + 1);
    sites_c = 1:1:(2*N);
    
    % 创建热图对比
    figure('Name', '各模型演化热图对比', 'Position', [100 100 1600 1200]);
    set(gcf,'Color','w','PaperPositionMode','auto');
    fontName = 'Helvetica';
    fontSize = 12;
    
    % 模型a的演化热图
    subplot(3, 2, 1);
    imagesc(sites_a_b, time_points, prob_a);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型a概率密度演化 (E_{scale}=%.3f)', E_scale_a),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    
    % 模型b的演化热图
    subplot(3, 2, 2);
    imagesc(sites_a_b, time_points, prob_b);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型b概率密度演化 (E_{scale}=%.3f)', E_scale_b),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    
    % 模型c（从N开始）的演化热图
    subplot(3, 2, 3);
    imagesc(sites_c, time_points, prob_c_N);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（从第N个格点开始）(E_{scale}=%.3f)', E_scale_c),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    
    % 模型c（从N+1开始）的演化热图
    subplot(3, 2, 4);
    imagesc(sites_c, time_points, prob_c_N1);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（从第N+1个格点开始）(E_{scale}=%.3f)', E_scale_c),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    
    % 模型c（从N和N+1各0.5开始）的演化热图
    subplot(3, 2, 5);
    imagesc(sites_c, time_points, prob_c_NN1);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（N和N+1各0.5）(E_{scale}=%.3f)', E_scale_c),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    
    % 添加总标题
    sgtitle(sprintf('各模型演化热图对比 (无量纲时间Et=%.1f)', total_time), 'FontSize', 16, 'FontWeight', 'bold');
    
    print('-dtiff','-r600','Models_Evolution.tiff');

    % 单独图像输出（保持原代码结构）
    figure('Name', 'Model-a', 'Position', [100 100 1200 800]);
    set(gcf,'Color','w','PaperPositionMode','auto');
    fontName = 'Helvetica';
    fontSize = 14;
    imagesc(sites_a_b, time_points, prob_a);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型a概率密度演化 (Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    print('-dtiff','-r600','Model-a.tiff');
    
    figure('Name', 'Model-b', 'Position', [100 100 1200 800]);
    set(gcf,'Color','w','PaperPositionMode','auto');
    fontName = 'Helvetica';
    fontSize = 14;
    imagesc(sites_a_b, time_points, prob_b);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型b概率密度演化 (Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    print('-dtiff','-r600','Model-b.tiff');
    
    figure('Name', 'Model-c-4', 'Position',[100 100 1200 800]);
    set(gcf,'Color','w','PaperPositionMode','auto');
    fontName = 'Helvetica';
    fontSize = 14;
    imagesc(sites_c, time_points, prob_c_N);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（从第N个格点开始）(Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    print('-dtiff','-r600','Model-c-4.tiff');
    
    figure('Name', 'Model-c-5', 'Position', [100 100 1200 800]);
    set(gcf,'Color','w','PaperPositionMode','auto');
    fontName = 'Helvetica';
    fontSize = 14;
    imagesc(sites_c, time_points, prob_c_N1);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（从第N+1个格点开始）(Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    print('-dtiff','-r600','Model-c-5.tiff');
    
    figure('Name', 'Model-c-NN1', 'Position', [100 100 1200 800]);
    set(gcf,'Color','w','PaperPositionMode','auto');
    fontName = 'Helvetica';
    fontSize = 14;
    imagesc(sites_c, time_points, prob_c_NN1);
    set(gca,'FontName',fontName,'FontSize',fontSize);
    colorbar;
    ylabel('无量纲时间 Et','FontSize',fontSize,'FontName',fontName);
    xlabel('格点位置','FontSize',fontSize,'FontName',fontName);
    title(sprintf('模型c（N和N+1各0.5）(Et=%.1f)', total_time),'FontSize',fontSize,'FontName',fontName);
    axis xy;
    grid on;
    box on;
    clim([0, 1]);
    print('-dtiff','-r600','Model-c-NN1.tiff');
end
% 模型a的函数声明（实际代码应包含完整实现）
function [H, eigenvalues, eigenvectors] = construct_hamiltonian_a(N, t, delta, gamma)
    % 构造模型a的哈密顿量并求解本征值和本征向量
    % 模型b与模型a的区别：后N个格点的增益耗散强度分布与模型a完全相反
    % N: 格点参数，总格点数为2N+1
    % t, delta: 跃迁强度参数
    % gamma: 增益耗散强度参数
    % H: 构造的哈密顿量
    % eigenvalues: 本征值
    % eigenvectors: 本征向量
    
    % 总格点数
    total_sites = 2*N + 1;
    
    % 初始化稀疏矩阵
    H = sparse(total_sites, total_sites);
    
    % 设置增益耗散强度（对角元素）
    for i = 1:total_sites
        if i == N + 1
            % 第N+1个格点的增益耗散强度为0（与模型a相同）
            H(i, i) = 0;
        else
            if mod(N, 2) == 1  % N为奇数
                % 模型a中后N个格点的增益耗散强度与模型b相反
                if i <= N  % 前N个格点，与模型a相同
                    if mod(i, 2) == 1  % 奇数位置
                        H(i, i) = -1i * gamma;
                    else  % 偶数位置
                        H(i, i) = 1i * gamma;
                    end
                else  % 后N个格点，与模型a相反
                    if mod(i, 2) == 1  % 奇数位置（模型a为-1i*gamma）
                        H(i, i) = 1i * gamma;
                    else  % 偶数位置（模型a为1i*gamma）
                        H(i, i) = -1i * gamma;
                    end
                end
            else  % N为偶数
                if i <= N  % 前N个格点，与模型a相同
                    if mod(i, 2) == 1  % 奇数位置
                        H(i, i) = 1i * gamma;
                    else  % 偶数位置
                        H(i, i) = -1i * gamma;
                    end
                else  % 后N个格点，与模型a相反
                    if mod(i, 2) == 1  % 奇数位置（模型a为-1i*gamma）
                        H(i, i) = 1i * gamma;
                    else  % 偶数位置（模型a为1i*gamma）
                        H(i, i) = -1i * gamma;
                    end
                end
            end
        end
    end
    
    % 设置跃迁强度（非对角元素）- 与模型a完全相同
    if mod(N, 2) == 1  % N为奇数
        % 前N个格点的跃迁
        % 第一个和第二个格点之间的跃迁强度为t-delta
        H(1, 2) = t - delta;
        H(2, 1) = conj(t - delta);  % 哈密顿量应该是厄米的
        
        % 剩下的N-1个格点，组成(N-1)/2个晶胞
        num_cells = (N - 1) / 2;
        for n = 1:num_cells
            % 第2n个和第2n+1个格点之间的跃迁强度为t+delta
            site1 = 2*n;
            site2 = 2*n + 1;
            if site2 <= N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            
            % 第2n+1个和第2n+2个格点之间的跃迁强度为t-delta
            site1 = 2*n + 1;
            site2 = 2*n + 2;
            if site2 <= N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
        
        % 连接前N个格点和中间格点(N+1)
        H(N, N+1) = t - delta;
        H(N+1, N) = conj(t - delta);
        
        % 后N个格点的跃迁 (从N+1到2N+1)
        for n = 1:num_cells
            % 第(N+1+2n-2)个和第(N+1+2n-1)个格点之间的跃迁强度为t-delta
            site1 = N + 1 + 2*n - 2;
            site2 = N + 1 + 2*n - 1;
            if site2 <= total_sites
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
            
            % 第(N+1+2n-1)个和第(N+1+2n)个格点之间的跃迁强度为t+delta
            site1 = N + 1 + 2*n - 1;
            site2 = N + 1 + 2*n;
            if site2 <= total_sites
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
        end
        
        % 最后一个和倒数第二个格点之间的跃迁强度为t-delta
        H(total_sites - 1, total_sites) = t - delta;
        H(total_sites, total_sites - 1) = conj(t - delta);
        
    else  % N为偶数
        % 前N个格点的跃迁
        num_cells = N / 2;
        for n = 1:num_cells
            % 第2n-1个和第2n个格点之间的跃迁强度为t+delta
            site1 = 2*n - 1;
            site2 = 2*n;
            if site2 <= N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            
            % 第2n个和第2n+1个格点之间的跃迁强度为t-delta
            site1 = 2*n;
            site2 = 2*n + 1;
            if site2 <= N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
        
        % 连接前N个格点和中间格点(N+1)
        H(N, N+1) = t - delta;
        H(N+1, N) = conj(t - delta);
        
        % 后N个格点的跃迁 (从N+1到2N+1)
        for n = 1:num_cells
            % 第(N+1+2n-2)个和第(N+1+2n-1)个格点之间的跃迁强度为t-delta
            site1 = N + 1 + 2*n - 2;
            site2 = N + 1 + 2*n - 1;
            if site1 >= N + 1 && site2 <= total_sites
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
            
            % 第(N+1+2n-1)个和第(N+1+2n)个格点之间的跃迁强度为t+delta
            site1 = N + 1 + 2*n - 1;
            site2 = N + 1 + 2*n;
            if site2 <= total_sites
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
        end
    end
    
    % 求解本征值和本征向量
    if total_sites <= 1000  % 对于小矩阵，使用full和eig
        H_full = full(H);
        [eigenvectors, eigenvalues_matrix] = eig(H_full);
        eigenvalues = diag(eigenvalues_matrix);
    else  % 对于大矩阵，使用稀疏矩阵的eigs
        % 这里求解全部本征值，对于非常大的矩阵可以只求解部分
        [eigenvectors, eigenvalues_matrix] = eigs(H, total_sites);
        eigenvalues = diag(eigenvalues_matrix);
    end
end

function [evolution, prob_density] = simulate_evolution_a(N, t, delta, gamma, total_time, num_steps)
    % 模拟模型a中初始状态在第N+1个格点的时间演化
    % N: 格点参数
    % t, delta, gamma: 模型参数
    % total_time: 总演化时间
    % num_steps: 时间步数
    % evolution: 演化结果，每行是一个时间步的状态
    % prob_density: 概率密度演化，每行是一个时间步的概率密度
    
    % 构造模型b的哈密顿量
    [H, ~, ~] = construct_hamiltonian_a(N, t, delta, gamma);
    total_sites = 2*N + 1;
    
    % 时间步长
    dt = total_time / num_steps;
    
    % 初始状态：第N+1个格点为1，其余为0
    initial_state = zeros(total_sites, 1);
    initial_state(N + 1) = 1;
    
    % 初始化演化结果矩阵
    evolution = zeros(num_steps + 1, total_sites);
    evolution(1, :) = initial_state';
    
    % 初始化概率密度矩阵
    prob_density = zeros(num_steps + 1, total_sites);
    prob_density(1, :) = abs(initial_state').^2;
    
    % 当前状态
    current_state = initial_state;
    
    % 时间演化
    for step = 1:num_steps
        % 计算时间演化算符：exp(-i*H*dt)
        if total_sites <= 1000
            % 小矩阵使用直接指数化
            U = expm(-1i * full(H) * dt);
        else
            % 大稀疏矩阵使用指数化的近似方法
            U = speye(total_sites) - 1i * H * dt - 0.5 * (H^2) * (dt^2);
        end
        
        % 应用时间演化算符
        current_state = U * current_state;
        
        % 归一化
        current_state = current_state / norm(current_state);
        
        % 保存当前状态和概率密度
        evolution(step + 1, :) = current_state';
        prob_density(step + 1, :) = abs(current_state').^2;
    end
end

% 模型b的函数声明（实际代码应包含完整实现）
function [H, eigenvalues, eigenvectors] = construct_hamiltonian_b(N, t, delta, gamma)
    % 构造模型b的哈密顿量并求解本征值和本征向量
    % N: 格点参数，总格点数为2N+1
    % t, delta: 跃迁强度参数
    % gamma: 增益耗散强度参数
    % H: 构造的哈密顿量
    % eigenvalues: 本征值
    % eigenvectors: 本征向量
    
    % 总格点数
    total_sites = 2*N + 1;
    
    % 初始化稀疏矩阵
    H = sparse(total_sites, total_sites);
    
    % 设置增益耗散强度（对角元素）
    for i = 1:total_sites
        if i == N + 1
            % 第N+1个格点的增益耗散强度为0
            H(i, i) = 0;
        else
            if mod(N, 2) == 1  % N为奇数
                if mod(i, 2) == 1  % 奇数位置
                    H(i, i) = -1i * gamma;
                else  % 偶数位置
                    H(i, i) = 1i * gamma;
                end
            else  % N为偶数
                if i <= N  % 前N个格点
                    if mod(i, 2) == 1  % 奇数位置
                        H(i, i) = 1i * gamma;
                    else  % 偶数位置
                        H(i, i) = -1i * gamma;
                    end
                else  % 后N个格点（不包括第N+1个）
                    if mod(i, 2) == 1  % 奇数位置
                        H(i, i) = -1i * gamma;
                    else  % 偶数位置
                        H(i, i) = 1i * gamma;
                    end
                end
            end
        end
    end
    
    % 设置跃迁强度（非对角元素）
    if mod(N, 2) == 1  % N为奇数
        % 前N个格点的跃迁
        % 第一个和第二个格点之间的跃迁强度为t-delta
        H(1, 2) = t - delta;
        H(2, 1) = conj(t - delta);  % 哈密顿量应该是厄米的
        
        % 剩下的N-1个格点，组成(N-1)/2个晶胞
        num_cells = (N - 1) / 2;
        for n = 1:num_cells
            % 第2n个和第2n+1个格点之间的跃迁强度为t+delta
            site1 = 2*n;
            site2 = 2*n + 1;
            if site2 <= N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            
            % 第2n+1个和第2n+2个格点之间的跃迁强度为t-delta
            site1 = 2*n + 1;
            site2 = 2*n + 2;
            if site2 <= N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
        
        % 连接前N个格点和中间格点(N+1)
        H(N, N+1) = t - delta;
        H(N+1, N) = conj(t - delta);
        
        % 后N个格点的跃迁 (从N+1到2N+1)
        for n = 1:num_cells
            % 第(N+1+2n-2)个和第(N+1+2n-1)个格点之间的跃迁强度为t-delta
            site1 = N + 1 + 2*n - 2;
            site2 = N + 1 + 2*n - 1;
            if site2 <= total_sites
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
            
            % 第(N+1+2n-1)个和第(N+1+2n)个格点之间的跃迁强度为t+delta
            site1 = N + 1 + 2*n - 1;
            site2 = N + 1 + 2*n;
            if site2 <= total_sites
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
        end
        
        % 最后一个和倒数第二个格点之间的跃迁强度为t-delta
        H(total_sites - 1, total_sites) = t - delta;
        H(total_sites, total_sites - 1) = conj(t - delta);
        
    else  % N为偶数
        % 前N个格点的跃迁
        num_cells = N / 2;
        for n = 1:num_cells
            % 第2n-1个和第2n个格点之间的跃迁强度为t+delta
            site1 = 2*n - 1;
            site2 = 2*n;
            if site2 <= N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            
            % 第2n个和第2n+1个格点之间的跃迁强度为t-delta
            site1 = 2*n;
            site2 = 2*n + 1;
            if site2 <= N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
        
        % 连接前N个格点和中间格点(N+1)
        H(N, N+1) = t - delta;
        H(N+1, N) = conj(t - delta);
        
        % 后N个格点的跃迁 (从N+1到2N+1)
        for n = 1:num_cells
            % 第(N+1+2n-2)个和第(N+1+2n-1)个格点之间的跃迁强度为t-delta
            site1 = N + 1 + 2*n - 2;
            site2 = N + 1 + 2*n - 1;
            if site1 >= N + 1 && site2 <= total_sites
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
            
            % 第(N+1+2n-1)个和第(N+1+2n)个格点之间的跃迁强度为t+delta
            site1 = N + 1 + 2*n - 1;
            site2 = N + 1 + 2*n;
            if site2 <= total_sites
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
        end
    end
    
    % 求解本征值和本征向量
    if total_sites <= 100  % 对于小矩阵，使用full和eig
        H_full = full(H);
        [eigenvectors, eigenvalues_matrix] = eig(H_full);
        eigenvalues = diag(eigenvalues_matrix);
    else  % 对于大矩阵，使用稀疏矩阵的eigs
        % 这里求解全部本征值，对于非常大的矩阵可以只求解部分
        [eigenvectors, eigenvalues_matrix] = eigs(H, total_sites);
        eigenvalues = diag(eigenvalues_matrix);
    end
end

function [evolution, prob_density] = simulate_evolution_b(N, t, delta, gamma, total_time, num_steps)
    % 模拟Model-b初始状态在第N+1个格点的时间演化
    % N: 格点参数
    % t, delta, gamma: 模型参数
    % total_time: 总演化时间
    % num_steps: 时间步数
    % evolution: 演化结果，每行是一个时间步的状态
    % prob_density: 概率密度演化，每行是一个时间步的概率密度
    
    % 构造哈密顿量
    [H, ~, ~] = construct_hamiltonian_b(N, t, delta, gamma);
    total_sites = 2*N + 1;
    
    % 时间步长
    dt = total_time / num_steps;
    
    % 初始状态：第N+1个格点为1，其余为0
    initial_state = zeros(total_sites, 1);
    initial_state(N + 1) = 1;
    
    % 初始化演化结果矩阵
    evolution = zeros(num_steps + 1, total_sites);
    evolution(1, :) = initial_state';
    
    % 初始化概率密度矩阵
    prob_density = zeros(num_steps + 1, total_sites);
    prob_density(1, :) = abs(initial_state').^2;
    
    % 当前状态
    current_state = initial_state;
    
    % 时间演化
    for step = 1:num_steps
        % 计算时间演化算符：exp(-i*H*dt)
        % 对于大矩阵，使用更高效的方法
        if total_sites <= 1000
            % 小矩阵使用直接指数化
            U = expm(-1i * full(H) * dt);
        else
            % 大稀疏矩阵使用指数化的近似方法
            % 采用泰勒展开到二阶
            U = speye(total_sites) - 1i * H * dt - 0.5 * (H^2) * (dt^2);
        end
        
        % 应用时间演化算符
        current_state = U * current_state;
        
        % 归一化
        current_state = current_state / norm(current_state);
        
        % 保存当前状态和概率密度
        evolution(step + 1, :) = current_state';
        prob_density(step + 1, :) = abs(current_state').^2;
    end
end

% Model-c
function [H, eigenvalues, eigenvectors] = construct_hamiltonian_c(N, t, delta, gamma)
    % 构造模型c的哈密顿量并求解本征值和本征向量
    % 模型c: 共有2N个格点，第N个与第N+1个格点间跃迁强度固定为单位强度t
    % N: 格点参数，总格点数为2N
    % t, delta: 跃迁强度参数
    % gamma: 增益耗散强度参数
    % H: 构造的哈密顿量
    % eigenvalues: 本征值
    % eigenvectors: 本征向量
    
    % 总格点数
    total_sites = 2*N;
    
    % 初始化稀疏矩阵
    H = sparse(total_sites, total_sites);
    
    % 设置增益耗散强度（对角元素）
    for i = 1:total_sites
        if i <= N  % 前N个格点
            if mod(N, 2) == 1  % N为奇数
                if mod(i, 2) == 1  % 奇数位置
                    H(i, i) = 1i * gamma;
                else  % 偶数位置
                    H(i, i) = -1i * gamma;
                end
            else  % N为偶数
                if mod(i, 2) == 1  % 奇数位置
                    H(i, i) = -1i * gamma;
                else  % 偶数位置
                    H(i, i) = 1i * gamma;
                end
            end
        else  % 后N个格点 (N+1到2N)
            % 后N个格点无论N奇偶，增益耗散分布相同
            if mod(i, 2) == 1  % 奇数位置
                H(i, i) = 1i * gamma;
            else  % 偶数位置
                H(i, i) = -1i * gamma;
            end
        end
    end
    
    % 设置跃迁强度（非对角元素）
    % 第N个和第N+1个格点之间的跃迁强度固定为单位强度t
    H(N, N+1) = t;
    H(N+1, N) = conj(t);  % 保持厄米性
    
    if mod(N, 2) == 1  % N为奇数
        % 前N个格点的跃迁 (1 ≤ n ≤ N)
        % 每两个格点组成原胞，共(N-1)/2个原胞
        num_cells = (N - 1) / 2;
        for cell = 1:num_cells
            % 原胞内的跃迁 (2cell-1 与 2cell)
            site1 = 2*cell - 1;
            site2 = 2*cell;
            if site2 <= N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            
            % 原胞间的跃迁 (2cell 与 2cell+1)
            site1 = 2*cell;
            site2 = 2*cell + 1;
            if site2 <= N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
        
        % 后N个格点的跃迁 (N+1 ≤ n ≤ 2N)
        % 第2N-1个与第2N个格点间的跃迁强度为t-delta
        H(2*N-1, 2*N) = t - delta;
        H(2*N, 2*N-1) = conj(t - delta);
        
        % 每两个格点组成原胞，共(N-1)/2个原胞
        for cell = 1:num_cells
            % 原胞内的跃迁
            site1 = N + 1 + 2*(cell - 1);
            site2 = N + 1 + 2*(cell - 1) + 1;
            if site2 <= 2*N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            
            % 原胞间的跃迁
            site1 = N + 1 + 2*(cell - 1) + 1;
            site2 = N + 1 + 2*(cell - 1) + 2;
            if site2 <= 2*N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
        
    else  % N为偶数
        % 前N个格点的跃迁 (1 ≤ n ≤ N)
        % 第1个与第2个格点间的跃迁强度为t-delta
        H(1, 2) = t - delta;
        H(2, 1) = conj(t - delta);
        
        % 第N-1个与第N个格点间的跃迁强度为t-delta
        H(N-1, N) = t - delta;
        H(N, N-1) = conj(t - delta);
        
        % 中间格点组成原胞，共(N-2)/2个原胞
        if N >= 2
            num_cells = (N - 2) / 2;
            for cell = 1:num_cells
                % 原胞内的跃迁
                site1 = 2 + 2*(cell - 1);
                site2 = 2 + 2*(cell - 1) + 1;
                if site2 <= N
                    H(site1, site2) = t + delta;
                    H(site2, site1) = conj(t + delta);
                end
                
                % 原胞间的跃迁
                site1 = 2 + 2*(cell - 1) + 1;
                site2 = 2 + 2*(cell - 1) + 2;
                if site2 <= N
                    H(site1, site2) = t - delta;
                    H(site2, site1) = conj(t - delta);
                end
            end
        end
        
        % 后N个格点的跃迁 (N+1 ≤ n ≤ 2N)
        % 每两个格点组成原胞，共N/2个原胞
        num_cells = N / 2;
        for cell = 1:num_cells
            % 原胞内的跃迁
            site1 = N + 1 + 2*(cell - 1);
            site2 = N + 1 + 2*(cell - 1) + 1;
            if site2 <= 2*N
                H(site1, site2) = t + delta;
                H(site2, site1) = conj(t + delta);
            end
            
            % 原胞间的跃迁
            site1 = N + 1 + 2*(cell - 1) + 1;
            site2 = N + 1 + 2*(cell - 1) + 2;
            if site2 <= 2*N
                H(site1, site2) = t - delta;
                H(site2, site1) = conj(t - delta);
            end
        end
    end
    
    % 求解本征值和本征向量
    if total_sites <= 1000  % 对于小矩阵，使用full和eig
        H_full = full(H);
        [eigenvectors, eigenvalues_matrix] = eig(H_full);
        eigenvalues = diag(eigenvalues_matrix);
    else  % 对于大矩阵，使用稀疏矩阵的eigs
        [eigenvectors, eigenvalues_matrix] = eigs(H, total_sites);
        eigenvalues = diag(eigenvalues_matrix);
    end
end

function [evolution_N, evolution_N1, prob_N, prob_N1, prob_NN1] = simulate_evolution_c(N, t, delta, gamma, total_time, num_steps)
     % 模拟模型c中初始状态在第N和N+1个格点的时间演化
    % 新增：第三种初态为第N和第N+1个格点上各为0.5
    % N: 格点参数
    % t, delta, gamma: 模型参数
    % total_time: 总演化时间
    % num_steps: 时间步数
    % evolution_N: 从第N个格点开始的演化结果
    % evolution_N1: 从第N+1个格点开始的演化结果
    % prob_N: 从第N个格点开始的概率密度演化
    % prob_N1: 从第N+1个格点开始的概率密度演化
    % prob_NN1: 从N和N+1各0.5开始的概率密度演化
    
    % 构造模型c的哈密顿量
    [H, ~, ~] = construct_hamiltonian_c(N, t, delta, gamma);
    total_sites = 2*N;
    
    % 时间步长
    dt = total_time / num_steps;
    
    % 初始化状态：分别在第N、第N+1个格点为1，以及N和N+1各0.5
    initial_state_N = zeros(total_sites, 1);
    initial_state_N(N) = 1;
    
    initial_state_N1 = zeros(total_sites, 1);
    initial_state_N1(N+1) = 1;
    
    initial_state_NN1 = zeros(total_sites, 1);
    initial_state_NN1(N) = 0.5;
    initial_state_NN1(N+1) = 0.5;
    initial_state_NN1 = initial_state_NN1 / norm(initial_state_NN1); % 归一化
    
    % 初始化演化结果矩阵
    evolution_N = zeros(num_steps + 1, total_sites);
    evolution_N(1, :) = initial_state_N';
    
    evolution_N1 = zeros(num_steps + 1, total_sites);
    evolution_N1(1, :) = initial_state_N1';
    
    evolution_NN1 = zeros(num_steps + 1, total_sites);
    evolution_NN1(1, :) = initial_state_NN1';
    
    % 初始化概率密度矩阵
    prob_N = zeros(num_steps + 1, total_sites);
    prob_N(1, :) = abs(initial_state_N').^2;
    
    prob_N1 = zeros(num_steps + 1, total_sites);
    prob_N1(1, :) = abs(initial_state_N1').^2;
    
    prob_NN1 = zeros(num_steps + 1, total_sites);
    prob_NN1(1, :) = abs(initial_state_NN1').^2;
    
    % 当前状态
    current_state_N = initial_state_N;
    current_state_N1 = initial_state_N1;
    current_state_NN1 = initial_state_NN1;
    
    % 时间演化
    for step = 1:num_steps
        % 计算时间演化算符：exp(-i*H*dt)
        if total_sites <= 1000
            % 小矩阵使用直接指数化
            U = expm(-1i * full(H) * dt);
        else
            % 大稀疏矩阵使用指数化的近似方法
            U = speye(total_sites) - 1i * H * dt - 0.5 * (H^2) * (dt^2);
        end
        
        % 应用时间演化算符
        current_state_N = U * current_state_N;
        current_state_N1 = U * current_state_N1;
        current_state_NN1 = U * current_state_NN1;
        
        % 归一化
        current_state_N = current_state_N / norm(current_state_N);
        current_state_N1 = current_state_N1 / norm(current_state_N1);
        current_state_NN1 = current_state_NN1 / norm(current_state_NN1);
        
        % 保存当前状态和概率密度
        evolution_N(step + 1, :) = current_state_N';
        evolution_N1(step + 1, :) = current_state_N1';
        evolution_NN1(step + 1, :) = current_state_NN1';
        
        prob_N(step + 1, :) = abs(current_state_N').^2;
        prob_N1(step + 1, :) = abs(current_state_N1').^2;
        prob_NN1(step + 1, :) = abs(current_state_NN1').^2;
    end
end

function w = expmv(t, A, v, tol)
    % 使用Krylov子空间方法计算 exp(t*A)*v
    % 这是一个简化版本，对于大型稀疏矩阵更高效
    
    if nargin < 4
        tol = 1e-10;
    end
    
    n = length(v);
    m = min(30, n); % Krylov子空间维数
    
    % 归一化初始向量
    beta = norm(v);
    v = v / beta;
    
    % Arnoldi过程
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
    
    % 计算 exp(t*H_m) * e_1
    e1 = zeros(m, 1);
    e1(1) = 1;
    expH_e1 = expm(t * H(1:m, 1:m)) * e1;
    
    % 重构解
    w = beta * V(:, 1:m) * expH_e1;
end
% 其他函数保持不变（construct_hamiltonian_a, construct_hamiltonian_b, construct_hamiltonian_c等）
% 由于代码长度限制，这里只显示修改的部分，其他函数保持原样