function [UU, V, A, Z, iter, obj] = Fdisp_qp(X, Y, lambda, numanchor, U,numclass)


    %% 初始化
    maxIter = 50;
    m = numanchor;  % 锚点数量，统一为k
    numview = length(X);
    numsample = size(Y, 1);  % 样本数，n

    % 使用传入的 U 替换 Z
    Z = U;
    disp('Initial size of Z:');
    disp(size(Z));  % 打印 Z 的初始尺寸

    % 初始化 A
    for i = 1:numview
        di = size(X{i}, 1); 
        A{i} = zeros(di, m);  % A 是 d_p x m 的大小
        disp(['Initial size of A{', num2str(i), '}: ', num2str(size(A{i}))]);  % 打印 A 的初始尺寸
    end

    alpha = ones(1, numview) / numview;

    flag = 1;
    iter = 0;

    %% 主要优化循环
    while flag
        iter = iter + 1;
        disp(['Iteration ', num2str(iter)]);

        %% 优化 A_p 矩阵
        for ia = 1:numview
            disp(['Size of X{', num2str(ia), '}: ', num2str(size(X{ia}))]);
            disp(['Size of Z: ', num2str(size(Z))]);
            
            % 确保 X{ia}' 和 Z 的维度匹配
            part1 = X{ia}' * Z;  % X{ia}' 是 [d_p, n]，Z 是 [n, m]
            disp(['Size of part1 (X{', num2str(ia), '}'' * Z): ', num2str(size(part1))]);

            % SVD 分解 part1
            [Unew, ~, Vnew] = svd(part1, 'econ');
            A{ia} = Unew * Vnew';  % 更新 A_p，保持 [d_p, m] 维度
            disp(['Updated size of A{', num2str(ia), '}: ', num2str(size(A{ia}))]);
        end

        %% 优化 Z（实际是优化传入的 U）
        C1 = 0;
        C2 = 0;
        for a = 1:numview
            C1 = C1 + alpha(a)^2;
            C2 = C2 + alpha(a)^2 * A{a}' * X{a}';  % X{a}' 是 [n, d_p]
            disp(['Size of A{', num2str(a), '}'' * X{', num2str(a), '}'' (C2 contribution): ', num2str(size(A{a}' * X{a}'))]);
        end
        C1 = C1 + lambda * ones(1, numsample);  % C1 是 [1, n]
        disp(['Size of C1: ', num2str(size(C1))]);
        disp(['Size of C2: ', num2str(size(C2))]);

        % 更新 Z
        for ii = 1:numsample
            idx = 1:numanchor;  % numanchor 对应特征数，即 Z 的列数
            ut = C2(idx, ii) ./ C1(ii);
            Z(ii, idx) = EProjSimplex_new(ut');  % 投影到单纯形，更新 Z 的第 ii 行
        end
        disp(['Updated size of Z: ', num2str(size(Z))]);

        %% 优化 alpha
        M = zeros(numview, 1);
        for iv = 1:numview
            disp(['Size of A{', num2str(iv), '}: ', num2str(size(A{iv}))]);
            disp(['Size of Z: ', num2str(size(Z))]);
            
            % 确保 A{iv} 和 Z 的维度匹配，不再转置 Z
            M(iv) = norm(X{iv}' - A{iv} * Z', 'fro')^2;  % A{iv} * Z 是 [d_p, m] * [n, m]
        end
        Mfra = M.^-1;
        Q = 1 / sum(Mfra);
        alpha = Q * Mfra;
        disp(['Updated alpha: ', num2str(alpha')]);

        %% 计算目标函数
        term1 = 0;
        term2 = 0;
        for iv = 1:numview
            term1 = term1 + alpha(iv)^2 * norm(X{iv}' - A{iv} * Z', 'fro')^2;  
        end
        term2 = lambda * norm(Z', 'fro')^2;
        obj(iter) = term1 + term2;
        disp(['Objective function value at iteration ', num2str(iter), ': ', num2str(obj(iter))]);

        %% 收敛条件检查
        if (iter > 1) && (abs((obj(iter-1) - obj(iter)) / obj(iter-1)) < 1e-3 || iter > maxIter || obj(iter) < 1e-10)
            [UU, ~, V] = svd(Z, 'econ');  % Z 是 [n, m]，不需要转置
            UU = UU(:, 1:numclass);
            disp('Converged. Final size of UU:');
            disp(size(UU));  % 打印最终 UU 的尺寸
            flag = 0;
        end
    end
end
