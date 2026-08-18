function [UU, V, A, Z, iter, obj] = F_qp(X, Y, lambda, numanchor, U, numclass)
    % F_qp: 优化过程的主函数
    % 输入:
    %   X: Cell array，包含每个视图的数据矩阵
    %   Y: 样本标签
    %   lambda: 正则化参数
    %   numanchor: 锚点数量
    %   U: 初始的 U 矩阵
    %   numclass: 亚型数量
    % 输出:
    %   UU: 优化后的左奇异矩阵
    %   V: 右奇异矩阵
    %   A: 视图的映射矩阵
    %   Z: 表示矩阵
    %   iter: 迭代次数
    %   obj: 目标函数值

    %% 初始化
    maxIter = 50;
    m = numanchor;  % 锚点数，与 demo_qp 一致
    numview = length(X);
    numsample = size(Y, 1);  % 样本数，n

    % 使用传入的 U 替换 Z
    Z = U;

    % 初始化 A
    for i = 1:numview
        di = size(X{i}, 1); 
        A{i} = zeros(di, m);  % A 是 d_p x m 的大小
    end

    alpha = ones(1, numview) / numview;

    flag = 1;
    iter = 0;

    %% 主要优化循环
    while flag
        iter = iter + 1;

        %% 优化 A_p 矩阵
        for ia = 1:numview
            part1 = X{ia}' * Z;  % X{ia}' 是 [d_p, n]，Z 是 [n, m]
            [Unew, ~, Vnew] = svd(part1, 'econ');
            A{ia} = Unew * Vnew';  % 更新 A_p，保持 [d_p, m] 维度
        end

        %% 优化 Z
        C1 = 0;
        C2 = 0;
        for a = 1:numview
            C1 = C1 + alpha(a)^2;
            C2 = C2 + alpha(a)^2 * A{a}' * X{a}';  % X{a}' 是 [n, d_p]
        end
        C1 = C1 + lambda * ones(1, numsample);  % 确保 C1 是 [1, n] 大小

        % 更新 Z
        for ii = 1:numsample
            idx = 1:numanchor;  % numanchor 对应特征数，即 Z 的列数
            ut = C2(idx, ii) ./ C1(ii);
            Z(ii, idx) = EProjSimplex_new(ut');  % 投影到单纯形，更新 Z 的第 ii 行
        end

        %% 优化 alpha
        M = zeros(numview, 1);
        for iv = 1:numview
            M(iv) = norm(X{iv}' - A{iv} * Z', 'fro')^2;  % A{iv} * Z 是 [d_p, m] * [n, m]
        end
        Mfra = M.^-1;
        Q = 1 / sum(Mfra);
        alpha = Q * Mfra;

        %% 计算目标函数
        term1 = 0;
        term2 = 0;
        for iv = 1:numview
            term1 = term1 + alpha(iv)^2 * norm(X{iv}' - A{iv} * Z', 'fro')^2;  
        end
        term2 = lambda * norm(Z', 'fro')^2;
        obj(iter) = term1 + term2;

        %% 收敛条件检查
        if (iter > 1) && (abs((obj(iter-1) - obj(iter)) / obj(iter-1)) < 1e-3 || iter > maxIter || obj(iter) < 1e-10)
            [UU, ~, V] = svd(Z, 'econ');  % Z 是 [n, m]，不需要转置
            UU = UU(:, 1:numclass);  % 提取前 numclass 列
            flag = 0;
        end
    end
end
