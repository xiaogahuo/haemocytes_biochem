function CH_index = custom_calinski_harabasz(X_data, pred_label)
    % 获取样本数和聚类数
    [N, d] = size(X_data);
    unique_labels = unique(pred_label);
    k = length(unique_labels);

    % 计算整体数据的质心
    overall_mean = mean(X_data, 1);

    % 计算簇间和簇内的方差
    Sb = 0;  % 簇间方差
    Sw = 0;  % 簇内方差
    for i = 1:k
        cluster_points = X_data(pred_label == unique_labels(i), :);
        cluster_size = size(cluster_points, 1);
        cluster_mean = mean(cluster_points, 1);

        % 簇间离散度
        Sb = Sb + cluster_size * norm(cluster_mean - overall_mean)^2;

        % 簇内离散度
        Sw = Sw + sum(vecnorm(cluster_points - cluster_mean, 2, 2).^2);
    end

    % 计算 Calinski-Harabasz 指数
    CH_index = (Sb / (k - 1)) / (Sw / (N - k));
end
