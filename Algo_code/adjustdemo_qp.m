function [UU, W, Z, Zi, iter, obj] = adjustdemo_qp(X, Y, A, AA, numanchor, beta)
    % m      : the number of anchor. the size of Z is m*n.
    % beta   : the hyper-parameter of regularization term.
    % X      : n*di

    %% Initialize
    maxIter = 50; % the number of iterations
    m = numanchor;  % Number of anchors
    numview = length(X);  % Number of views
    numsample = size(Y, 1);  % Number of sample
    
    Z = zeros(m, numsample);  % Initialize Z (m * n)
    
    for i = 1:numview
        di = size(X{i}, 1);  % Dimension of view i
        W{i} = eye(m);  % Initialize W as an identity matrix
        X{i} = X{i}';  % Turn X{i} into d * n
        Zi{i} = zeros(m, numsample);  % Initialize Zi
    end
    Z(:, 1:m) = eye(m);  % Initialize Z as identity for the first m samples
    
    flag = 1;
    iter = 0;
    
    %% Main loop
    while flag
        iter = iter + 1;
        
        %% Optimize Z
        for iv = 1:numview
            H = W{iv} * AA{iv} * W{iv};  % m x m
            H = (H + H') / 2;  % Symmetrize H
            options = optimset('Algorithm', 'interior-point-convex', 'Display', 'off');
            Zp = zeros(m, numsample);  % Initialize Zp
            
            parfor j = 1:numsample
                ff = -X{iv}(:, j)' * A{iv} * W{iv};  % d * 1' * (d * m) * (m * m) = 1 * m
                Zp(:, j) = quadprog(H, ff', [], [], ones(1, m), 1, zeros(m, 1), ones(m, 1), [], options);
            end
            Zi{iv} = Zp;  % Save the updated Zp for view iv
            zz{iv} = Zp * Zp';  % Update zz{iv}
        end
        
        %% Optimize W_i
        options = optimset('Algorithm', 'interior-point-convex', 'Display', 'off');
        parfor iv = 1:numview
            R = zz{iv} .* AA{iv} + beta * AA{iv};  % Element-wise multiplication zz and AA
            fw = zeros(m, 1);
            G = Zi{iv} * X{iv}' * A{iv};  % G should be m x m
            for j = 1:m
                fw(j) = -G(j, j);  % Diagonal of G
            end
            ww{iv} = quadprog(R, fw', [], [], ones(1, m), m, zeros(m, 1), m * ones(m, 1), [], options);
        end
        
        % Update W{iv}
        for iv = 1:numview
            for j = 1:m
                W{iv}(j, j) = ww{iv}(j);
            end
        end
        
        %% Compute objective
        term1 = 0;
        term2 = 0;
        for iv = 1:numview
            term1 = term1 + norm(X{iv} - A{iv} * W{iv} * Zi{iv}, 'fro') ^ 2;
            term2 = term2 + ww{iv}' * AA{iv} * ww{iv};
        end
        obj(iter) = term1 + beta * term2;
        
        %% Check convergence
        if (iter > 1) && (abs((obj(iter - 1) - obj(iter)) / obj(iter - 1)) < 1e-3 || iter > maxIter || obj(iter) < 1e-10)
            Z = [];
            for iv = 1:numview
                Z = cat(1, Z, 1/sqrt(numview) * Zi{iv});
            end
            [UU, ~, ~] = mySVD(Z', numanchor);  % 调整为numanchor，维度由numanchor控制
            fprintf('The dimensions of UU matrix are: [%d, %d]\n', size(UU, 1), size(UU, 2));
            flag = 0;
        end
    end
end
