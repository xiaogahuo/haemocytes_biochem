function [X1, ind] = findindex(data, index)
%FINDINDEX Summary of this function goes here
%   Detailed explanation goes here
[numofview,~] = size(data);
[~,numofsample] = size(data{1});

X1 = cell(numofview,1);

ind = zeros(numofsample,numofview);
for i=1:numofview
    [d,~]=size(data{i});
    ind(index{i}, i) = 1;
    origin = data{i};
    origin(isnan(origin)) = 0;
    X1{i} = NormalizeData(origin);
 fprintf('View %d: Size of data{%d}: [%d, %d], Size of index{%d}: [%d], Non-zero elements in ind: %d\n', ...
    i, i, size(data{i}, 1), size(data{i}, 2), i, length(index{i}), nnz(ind(:, i)));

end
fprintf('Size of ind after update for view %d: [%d, %d]\n', i, size(ind, 1), size(ind, 2));


end

