clear;
clc;
warning off;
addpath(genpath('./'));

%% ==================== 参数设置 ====================
subtype_k_range   = [2,3,4,5,6,7];              % 亚型数量范围
anchor_multiplier = [1,2,3,5];                  % 锚点数 = subtype_k * multiplier
beta_range        = [2^(-2), 2^1, 2^4, 2^7];    % β 参数范围
lambda_grid       = [0.001, 0.1, 1, 10, 100];   % λ 参数范围

root_path = "/data4/myuan/"
path = root_path +"Hae_Bio/"
data_path = path + "dataset/"
multiview_result_path= path + "multiview_cluster_0/"

ds = {'Hae_Bio_0.0'};
%% ==================== 主循环 ====================
for dsi = 1:length(ds)
    dataName = ds{dsi};
    disp(['Processing dataset: ', dataName]);
    load(strcat(data_path, dataName));  % 加载 X, Y, eid

    numview = length(X);
    %% ---------- 遍历 subtype_k -----------
    for subtype_k = subtype_k_range
        fprintf('\n========== Running subtype_k = %d ==========\n', subtype_k);

        % 锚点范围
        anchor_range = subtype_k * anchor_multiplier;
        % 网格搜索结果
        csv_metrics = fullfile(multiview_result_path, sprintf('%s_Subtype_%d_clustering_results.csv', dataName, subtype_k));
        fid = fopen(csv_metrics, 'w');
        header = {'Anchor','Beta','Lambda',...
                  'Initial Silhouette','Optimized Silhouette',...
                  'Initial CH','Optimized CH','Initial DB','Optimized DB','Time'};
        fprintf(fid, '%s,', header{1:end-1});
        fprintf(fid, '%s\n', header{end});
        fclose(fid);

        %% ---------- 遍历 beta, anchor, lambda ----------
        for beta = beta_range
            for anchor = anchor_range

                % 初始化锚点
                X_norm = cell(size(X));
                A      = cell(size(X));
                AA     = cell(size(X));

                for v = 1:numview
                    rand('twister', 12);
                    X_norm{v} = mapstd(X{v}', 0, 1)';
                    [~, A{v}] = litekmeans(X_norm{v}, anchor, 'maxiter', 100, 'replicates', 3);
                    AA{v} = A{v} * A{v}';
                    A{v}  = A{v}';
                end

                % === Step 1：adjustdemo_qp → 初始 U
                tic;
                [U, W, Z, iter, obj] = adjustdemo_qp(X_norm, Y, A, AA, anchor, beta);
                timer1 = toc;

                [pred_label_init, ~] = myNMIACCwithmean(U, subtype_k);
                silhouette_init = mean(silhouette(U, pred_label_init));
                ch_init = custom_calinski_harabasz(U, pred_label_init);
                db_init = custom_davies_bouldin(U, pred_label_init);

                % === Step 2：遍历 lambda → F_qp 优化
                for lambda = lambda_grid

                    tic;
                    [UU_new, W_new, A_new, Z_new, iter_new, obj_new] = ...
                        F_qp(X_norm, Y, lambda, anchor, U, subtype_k);
                    timer2 = toc;

                    [pred_label_new, ~] = myNMIACCwithmean(UU_new, subtype_k);
                    silhouette_new = mean(silhouette(UU_new, pred_label_new));
                    ch_new = custom_calinski_harabasz(UU_new, pred_label_new);
                    db_new = custom_davies_bouldin(UU_new, pred_label_new);

                    total_time = timer1 + timer2;

                    % ---- 写入 grid-search CSV（文件2） ----
                    fid = fopen(csv_metrics, 'a');
                    fprintf(fid, '%d, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f\n',...
                        anchor, beta, lambda, ...
                        silhouette_init, silhouette_new, ...
                        ch_init, ch_new, db_init, db_new, total_time);
                    fclose(fid);

                    fprintf('Saved metrics: k=%d, a=%d, b=%.4f, l=%.4f\n', ...
                            subtype_k, anchor, beta, lambda);

                    % ---- 生成eid cluster_num CSV（每个参数组合一份） ----
                    csv_res = fullfile(multiview_result_path, ...
                        sprintf('%s_Subtype_%d_anchor_%d_beta_%.4f_lambda_%.4f_results.csv', ...
                            dataName, subtype_k, anchor, beta, lambda));

                    cluster_table = table(eid, pred_label_new, ...
                        'VariableNames', {'eid','final_cluster'});

                    writetable(cluster_table, csv_res);
                    fprintf('➡ 保存结果文件: %s\n', csv_res);

                end % lambda
            end % anchor
        end % beta
    end % subtype_k
end % dataset

disp('All parameter combinations complete! All results (metrics + file1-style CSVs) are saved.');
