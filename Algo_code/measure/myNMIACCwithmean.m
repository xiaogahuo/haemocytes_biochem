function [pred_label, std_val] = myNMIACCwithmean(U, numclass)

% 归一化
stream = RandStream.getGlobalStream;
reset(stream);
U_normalized = U ./ repmat(sqrt(sum(U.^2, 2)), 1, size(U, 2));  % 数据归一化
disp(['The dimensions of U_normalized matrix are: ', num2str(size(U_normalized))]);

maxIter = 50;

% 初始化结果存储
labels = zeros(size(U, 1), maxIter);  % 用来存储所有迭代的标签

for iter = 1:maxIter
    % 使用指定的亚型数量进行聚类，确保输入参数完整
    indx = litekmeans(U_normalized, numclass, 'distance', 'sqeuclidean', 'maxiter', 100, 'replicates', 3);  
    labels(:, iter) = indx(:);  % 保存每次迭代的标签
end

% 返回最终标签为最后一次迭代结果
pred_label = labels(:, end); 

% 计算标签的标准差
std_val = std(double(labels), 0, 2);  % 计算每个数据点标签在不同迭代中的标准差

end
