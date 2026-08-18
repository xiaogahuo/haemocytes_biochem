function DB_index = custom_davies_bouldin(X_data, pred_label)
    % 获取样本数和聚类数
    unique_labels = unique(pred_label);
    k = length(unique_labels);
    
    % 计算每个簇的质心
    centroids = zeros(k, size(X_data, 2));
    for i = 1:k
        cluster_points = X_data(pred_label == unique_labels(i), :);
        centroids(i, :) = mean(cluster_points, 1);  % 质心
    end
    
    % 计算每个簇的簇内距离（紧密度）
    sigma = zeros(k, 1);
    for i = 1:k
        cluster_points = X_data(pred_label == unique_labels(i), :);
        sigma(i) = mean(vecnorm(cluster_points - centroids(i, :), 2, 2));
    end
    
    % 计算簇间距离
    db_vals = zeros(k, 1);
    for i = 1:k
        max_ratio = -Inf;
        for j = 1:k
            if i ~= j
                d_ij = norm(centroids(i, :) - centroids(j, :));  % 簇中心之间的距离
                ratio = (sigma(i) + sigma(j)) / d_ij;  % 紧密度与簇间距离的比值
                max_ratio = max(max_ratio, ratio);
            end
        end
        db_vals(i) = max_ratio;  % 对于每个簇，选取最大比值
    end
    
    % 计算最终的 Davies-Bouldin 指数
    DB_index = mean(db_vals);
end
