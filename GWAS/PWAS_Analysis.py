import os
import numpy as np
import pandas as pd
import subprocess
import tempfile
import warnings
from scipy import stats
from scipy.stats import chi2
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
warnings.filterwarnings('ignore')

class PWASAnalyzer:
    def __init__(self, weights_dir, ldref_dir, output_dir, plink_bin='plink'):
        self.weights_dir = weights_dir
        self.ldref_dir = ldref_dir
        self.output_dir = output_dir
        self.plink_bin = plink_bin
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(os.path.join(output_dir, 'tmp'), exist_ok=True)

        self.pos_file = os.path.join(weights_dir, 'train_weights.pos.addN')
        self.pos_df = pd.read_csv(self.pos_file, sep='\t')
        print(f"[INFO] 加载权重位置文件: {self.pos_df.shape[0]} 个基因模型")

    def read_rdat_weights(self, rdat_path):
        try:
            import rpy2.robjects as ro
            from rpy2.robjects import pandas2ri
            pandas2ri.activate()
            ro.r(f'load("{rdat_path}")')
            wgt_obj = ro.r('wgt')
            if isinstance(wgt_obj, ro.vectors.ListVector):
                wgt_dict = {}
                for name in wgt_obj.names:
                    obj = wgt_obj[name]
                    if hasattr(obj, 'shape'):
                        wgt_dict[name] = np.array(obj)
                    else:
                        wgt_dict[name] = np.array(obj)
                if 'wgt' in wgt_dict:
                    wgt_data = wgt_dict['wgt']
                    if wgt_data.ndim == 2:
                        return wgt_data[:, 0] if wgt_data.shape[1] >= 1 else wgt_data.flatten()
                    return wgt_data.flatten()
                for key in ['top', 'bslmm', 'enet', 'lasso']:
                    if key in wgt_dict:
                        w = wgt_dict[key]
                        if w.ndim == 2:
                            return w[:, 0] if w.shape[1] >= 1 else w.flatten()
                        return w.flatten()
                if len(wgt_dict) > 0:
                    first_key = list(wgt_dict.keys())[0]
                    w = wgt_dict[first_key]
                    return w.flatten() if w.ndim > 1 else w
            return np.array(wgt_obj).flatten()
        except Exception as e:
            print(f"[WARN] rpy2 读取失败 ({e}), 尝试使用 Rscript 转换...")
            return self._read_rdat_via_rscript(rdat_path)

    def _read_rdat_via_rscript(self, rdat_path):
        tmp_csv = os.path.join(self.output_dir, 'tmp', 'wgt_tmp.csv')
        r_script = f'''
load("{rdat_path}")
if (is.list(wgt)) {{
    if ("wgt" %in% names(wgt)) {{
        w <- wgt[["wgt"]]
    }} else {{
        w <- wgt[[1]]
    }}
}} else {{
    w <- wgt
}}
if (is.matrix(w)) {{
    write.csv(data.frame(SNP=rownames(w), weight=w[,1], stringsAsFactors=FALSE), "{tmp_csv}", row.names=FALSE)
}} else {{
    write.csv(data.frame(SNP=names(w), weight=as.numeric(w), stringsAsFactors=FALSE), "{tmp_csv}", row.names=FALSE)
}}
'''
        r_file = os.path.join(self.output_dir, 'tmp', 'read_wgt.R')
        with open(r_file, 'w') as f:
            f.write(r_script)
        try:
            subprocess.run(['Rscript', r_file], check=True, capture_output=True, text=True)
            wgt_df = pd.read_csv(tmp_csv)
            return wgt_df['weight'].values
        except Exception as e:
            print(f"[ERROR] Rscript 转换也失败: {e}")
            return None

    def read_rdat_snp_info(self, rdat_path):
        return self._read_snp_info_via_rscript(rdat_path)

    def _read_snp_info_via_rscript(self, rdat_path):
        tmp_csv = os.path.join(self.output_dir, 'tmp', 'snp_info_tmp.csv')
        r_script = f'''
load("{rdat_path}")
if (is.list(wgt)) {{
    if ("wgt" %in% names(wgt)) {{
        w <- wgt[["wgt"]]
    }} else {{
        w <- wgt[[1]]
    }}
}} else {{
    w <- wgt
}}
if (is.matrix(w) && ncol(w) >= 4) {{
    df <- data.frame(SNP=w[,1], A1=w[,2], A2=w[,3], weight=w[,4], stringsAsFactors=FALSE)
}} else if (is.matrix(w) && ncol(w) >= 2) {{
    df <- data.frame(SNP=rownames(w), A1="", A2="", weight=w[,1], stringsAsFactors=FALSE)
}} else {{
    df <- data.frame(SNP=names(w), A1="", A2="", weight=as.numeric(w), stringsAsFactors=FALSE)
}}
write.csv(df, "{tmp_csv}", row.names=FALSE)
'''
        r_file = os.path.join(self.output_dir, 'tmp', 'read_snp.R')
        with open(r_file, 'w') as f:
            f.write(r_script)
        try:
            subprocess.run(['Rscript', r_file], check=True, capture_output=True, text=True)
            return pd.read_csv(tmp_csv)
        except Exception as e:
            print(f"[ERROR] 读取 SNP 信息失败: {e}")
            return None

    def compute_ld_matrix(self, chr_num, snp_list, window_bp=500000):
        tmp_dir = os.path.join(self.output_dir, 'tmp')
        bed_prefix = os.path.join(self.ldref_dir, f'1000G.EUR.{chr_num}')

        snp_file = os.path.join(tmp_dir, f'ld_snps_chr{chr_num}.txt')
        with open(snp_file, 'w') as f:
            for s in snp_list:
                f.write(s + '\n')

        extract_prefix = os.path.join(tmp_dir, f'ld_extract_chr{chr_num}')
        cmd = (
            f'{self.plink_bin} --bfile {bed_prefix} '
            f'--extract {snp_file} '
            f'--make-bed --out {extract_prefix} '
            f'--allow-no-sex'
        )
        try:
            subprocess.run(cmd.split(), check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            print(f"[WARN] PLINK extract 失败 chr{chr_num}: {e.stderr[:200]}")
            return None

        ld_prefix = os.path.join(tmp_dir, f'ld_matrix_chr{chr_num}')
        cmd_ld = (
            f'{self.plink_bin} --bfile {extract_prefix} '
            f'--r square --out {ld_prefix} '
            f'--allow-no-sex'
        )
        try:
            subprocess.run(cmd_ld.split(), check=True, capture_output=True, text=True)
            ld_file = ld_prefix + '.ld'
            if os.path.exists(ld_file):
                ld_mat = np.loadtxt(ld_file)
                if ld_mat.ndim == 1:
                    n_snps = int(np.sqrt(len(ld_mat)))
                    ld_mat = ld_mat[:n_snps * n_snps].reshape(n_snps, n_snps)
                return ld_mat
        except subprocess.CalledProcessError as e:
            print(f"[WARN] PLINK LD 计算失败 chr{chr_num}: {e.stderr[:200]}")
        except Exception as e:
            print(f"[WARN] LD 矩阵读取失败 chr{chr_num}: {e}")
        return None

    def compute_ld_from_bfile(self, chr_num, snp_list):
        tmp_dir = os.path.join(self.output_dir, 'tmp')
        bed_prefix = os.path.join(self.ldref_dir, f'1000G.EUR.{chr_num}')

        snp_file = os.path.join(tmp_dir, f'ld_snps_chr{chr_num}.txt')
        with open(snp_file, 'w') as f:
            for s in snp_list:
                f.write(s + '\n')

        extract_prefix = os.path.join(tmp_dir, f'ld_extract2_chr{chr_num}')
        cmd = (
            f'{self.plink_bin} --bfile {bed_prefix} '
            f'--extract {snp_file} '
            f'--make-bed --out {extract_prefix} '
            f'--allow-no-sex'
        )
        try:
            subprocess.run(cmd.split(), check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError:
            return None, None

        bim_file = extract_prefix + '.bim'
        if not os.path.exists(bim_file):
            return None, None
        bim_df = pd.read_csv(bim_file, sep='\t', header=None,
                             names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])

        ld_out = os.path.join(tmp_dir, f'ld_r_chr{chr_num}')
        cmd_ld = (
            f'{self.plink_bin} --bfile {extract_prefix} '
            f'--r --out {ld_out} '
            f'--allow-no-sex'
        )
        try:
            subprocess.run(cmd_ld.split(), check=True, capture_output=True, text=True)
            ld_file = ld_out + '.ld'
            if os.path.exists(ld_file):
                ld_df = pd.read_csv(ld_file, sep='\t', header=None,
                                    names=['SNP_A', 'SNP_B', 'R'])
                n_snps = len(bim_df)
                ld_mat = np.eye(n_snps)
                snp_names = bim_df['SNP'].values
                snp_to_idx = {s: i for i, s in enumerate(snp_names)}
                for _, row in ld_df.iterrows():
                    i = snp_to_idx.get(row['SNP_A'])
                    j = snp_to_idx.get(row['SNP_B'])
                    if i is not None and j is not None:
                        ld_mat[i, j] = row['R']
                        ld_mat[j, i] = row['R']
                return ld_mat, bim_df
        except Exception as e:
            print(f"[WARN] LD 计算失败 chr{chr_num}: {e}")
        return None, None

    def run_pwas(self, gwas_sumstats_path, snp_col='SNP', chr_col='CHR',
                 bp_col='BP', a1_col='A1', a2_col='A2',
                 beta_col='Beta', se_col='SE', p_col='P',
                 use_or=False, or_col='OR'):
        print("\n" + "=" * 60)
        print("开始 PWAS 分析")
        print("=" * 60)

        gwas_df = pd.read_csv(gwas_sumstats_path, sep=r'\s+|,', engine='python')
        print(f"[INFO] 加载 GWAS 汇总统计: {gwas_df.shape[0]} 个 SNP")

        gwas_df['Z'] = gwas_df[beta_col].astype(float) / gwas_df[se_col].astype(float)

        gwas_df[snp_col] = gwas_df[snp_col].astype(str)
        gwas_df[chr_col] = gwas_df[chr_col].astype(int)
        gwas_dict = {}
        for _, row in gwas_df.iterrows():
            gwas_dict[row[snp_col]] = {
                'Z': row['Z'], 'A1': str(row[a1_col]).upper(),
                'A2': str(row[a2_col]).upper(), 'CHR': row[chr_col],
                'BP': row[bp_col], 'Beta': row[beta_col], 'SE': row[se_col]
            }
        print(f"[INFO] GWAS Z-score 范围: [{gwas_df['Z'].min():.3f}, {gwas_df['Z'].max():.3f}]")

        results = []
        total_genes = self.pos_df.shape[0]

        for idx, row in self.pos_df.iterrows():
            gene_id = row['ID']
            chr_num = int(row['CHR'])
            wgt_path = os.path.join(self.weights_dir, row['WGT'])
            n_train = row.get('N', 0)

            if (idx + 1) % 100 == 0 or idx == 0:
                print(f"[INFO] 处理进度: {idx + 1}/{total_genes} ({gene_id})")

            snp_info = self.read_rdat_snp_info(wgt_path)
            if snp_info is None:
                weights = self.read_rdat_weights(wgt_path)
                if weights is None:
                    continue
                results.append({
                    'ID': gene_id, 'CHR': chr_num,
                    'P0': row['P0'], 'P1': row['P1'], 'N_train': n_train,
                    'n_snps': len(weights), 'Z': np.nan, 'P': np.nan,
                    'status': 'no_snp_info'
                })
                continue

            if isinstance(snp_info, pd.DataFrame) and 'SNP' in snp_info.columns:
                wgt_snps = snp_info['SNP'].astype(str).values
                wgt_a1 = snp_info['A1'].astype(str).str.upper().values if 'A1' in snp_info.columns else None
                wgt_a2 = snp_info['A2'].astype(str).str.upper().values if 'A2' in snp_info.columns else None
                wgt_vals = snp_info['weight'].astype(float).values if 'weight' in snp_info.columns else None
            else:
                wgt_snps = snp_info[:, 0].astype(str) if snp_info.shape[1] > 0 else np.array([])
                wgt_a1 = snp_info[:, 1].astype(str).str.upper() if snp_info.shape[1] > 1 else None
                wgt_a2 = snp_info[:, 2].astype(str).str.upper() if snp_info.shape[1] > 2 else None
                wgt_vals = snp_info[:, 3].astype(float) if snp_info.shape[1] > 3 else None

            if wgt_vals is None:
                weights = self.read_rdat_weights(wgt_path)
                if weights is None or len(weights) != len(wgt_snps):
                    continue
                wgt_vals = weights

            matched_indices = []
            matched_weights = []
            matched_gwas_z = []

            for i, snp in enumerate(wgt_snps):
                if snp in gwas_dict:
                    gwas_info = gwas_dict[snp]
                    w = wgt_vals[i]
                    z = gwas_info['Z']

                    if wgt_a1 is not None and wgt_a2 is not None:
                        if (wgt_a1[i] == gwas_info['A2'] and wgt_a2[i] == gwas_info['A1']):
                            w = -w

                    matched_indices.append(i)
                    matched_weights.append(w)
                    matched_gwas_z.append(z)

            if len(matched_weights) < 2:
                results.append({
                    'ID': gene_id, 'CHR': chr_num,
                    'P0': row['P0'], 'P1': row['P1'], 'N_train': n_train,
                    'n_snps': len(wgt_snps), 'n_matched': len(matched_weights),
                    'Z': np.nan, 'P': np.nan, 'status': 'too_few_snps'
                })
                continue

            w_vec = np.array(matched_weights)
            z_vec = np.array(matched_gwas_z)
            matched_snps = wgt_snps[matched_indices]

            ld_mat, bim_df = self.compute_ld_from_bfile(chr_num, matched_snps.tolist())

            if ld_mat is not None and ld_mat.shape[0] == len(w_vec):
                R = ld_mat
            else:
                n = len(w_vec)
                R = np.eye(n)

            R = np.clip(R, -1, 1)
            np.fill_diagonal(R, 1.0)

            try:
                eigenvalues = np.linalg.eigvalsh(R)
                if np.any(eigenvalues < -1e-10):
                    R = R + (abs(min(eigenvalues.min(), 0)) + 1e-6) * np.eye(len(R))

                numerator = w_vec @ z_vec
                variance = w_vec @ R @ w_vec
                if variance <= 0:
                    variance = w_vec @ w_vec

                z_pwas = numerator / np.sqrt(variance)
                p_pwas = 2 * stats.norm.sf(abs(z_pwas))

                results.append({
                    'ID': gene_id, 'CHR': chr_num,
                    'P0': row['P0'], 'P1': row['P1'], 'N_train': n_train,
                    'n_snps': len(wgt_snps), 'n_matched': len(matched_weights),
                    'Z': z_pwas, 'P': p_pwas,
                    'status': 'success'
                })
            except Exception as e:
                results.append({
                    'ID': gene_id, 'CHR': chr_num,
                    'P0': row['P0'], 'P1': row['P1'], 'N_train': n_train,
                    'n_snps': len(wgt_snps), 'n_matched': len(matched_weights),
                    'Z': np.nan, 'P': np.nan, 'status': f'error: {str(e)[:50]}'
                })

        result_df = pd.DataFrame(results)
        success_mask = result_df['status'] == 'success'
        if success_mask.sum() > 0:
            valid_p = result_df.loc[success_mask, 'P'].dropna()
            if len(valid_p) > 0:
                reject, fdr_pvals, _, _ = multipletests(valid_p, method='fdr_bh')
                result_df.loc[success_mask, 'P_FDR'] = fdr_pvals
                result_df.loc[success_mask, 'Significant_FDR0.05'] = reject

                bonf_pvals, _, _, _ = multipletests(valid_p, method='bonferroni')
                result_df.loc[success_mask, 'P_Bonferroni'] = bonf_pvals

        result_file = os.path.join(self.output_dir, 'PWAS_results.csv')
        result_df.to_csv(result_file, index=False)
        print(f"\n[INFO] PWAS 结果已保存至: {result_file}")
        print(f"[INFO] 总基因数: {len(results)}, 成功: {success_mask.sum()}")

        if success_mask.sum() > 0:
            sig_count = result_df.loc[success_mask, 'P'].dropna()
            print(f"[INFO] P < 0.05: {(sig_count < 0.05).sum()}")
            print(f"[INFO] FDR < 0.05: {result_df.loc[success_mask, 'Significant_FDR0.05'].sum()}")
            print(f"[INFO] P < 5e-6 (基因级别显著): {(sig_count < 5e-6).sum()}")

        return result_df

    def plot_manhattan(self, result_df, output_file=None, sig_level=0.05):
        if output_file is None:
            output_file = os.path.join(self.output_dir, 'PWAS_Manhattan.png')

        df = result_df.dropna(subset=['P', 'CHR']).copy()
        df = df[df['CHR'] != 0]
        df['CHR'] = df['CHR'].astype(int)
        df['-log10P'] = -np.log10(df['P'].clip(lower=1e-300))

        fig, ax = plt.subplots(figsize=(18, 6))
        colors = ['#1f77b4', '#ff7f0e']
        offsets = [0]
        max_log_p = 0

        for i, chr_num in enumerate(sorted(df['CHR'].unique())):
            chr_df = df[df['CHR'] == chr_num]
            if len(chr_df) == 0:
                continue
            x_pos = np.arange(len(chr_df)) + offsets[-1]
            ax.scatter(x_pos, chr_df['-log10P'].values,
                       c=colors[i % 2], s=8, alpha=0.7, edgecolors='none')
            max_log_p = max(max_log_p, chr_df['-log10P'].max())
            offsets.append(offsets[-1] + len(chr_df) + 5)

        sig_threshold = -np.log10(sig_level / df.shape[0]) if df.shape[0] > 0 else 0
        ax.axhline(y=sig_threshold, color='red', linestyle='--', linewidth=0.8,
                    label=f'Bonferroni (p={sig_level}/{df.shape[0]})')
        ax.axhline(y=-np.log10(0.05), color='blue', linestyle=':', linewidth=0.6,
                    label='Nominal (p=0.05)')

        tick_positions = []
        tick_labels = []
        offset_idx = 0
        for chr_num in sorted(df['CHR'].unique()):
            chr_df = df[df['CHR'] == chr_num]
            if len(chr_df) > 0:
                tick_positions.append(offset_idx + len(chr_df) / 2)
                tick_labels.append(str(chr_num))
                offset_idx += len(chr_df) + 5

        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, fontsize=8)
        ax.set_xlabel('Chromosome', fontsize=12)
        ax.set_ylabel('-log10(P)', fontsize=12)
        ax.set_title('PWAS Manhattan Plot', fontsize=14)
        ax.legend(fontsize=9)
        ax.set_ylim(0, min(max_log_p * 1.1, 50))

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] Manhattan plot 已保存: {output_file}")

    def plot_qq(self, result_df, output_file=None):
        if output_file is None:
            output_file = os.path.join(self.output_dir, 'PWAS_QQ.png')

        df = result_df.dropna(subset=['P']).copy()
        df = df[df['P'] > 0]
        if len(df) == 0:
            print("[WARN] 无有效 P 值用于 QQ plot")
            return

        observed = -np.log10(df['P'].sort_values().values)
        n = len(observed)
        expected = -np.log10(np.arange(1, n + 1) / (n + 1))

        fig, ax = plt.subplots(figsize=(7, 7))
        ax.scatter(expected, observed, s=10, c='#1f77b4', alpha=0.7, edgecolors='none')
        ax.plot([0, max(expected) + 0.5], [0, max(expected) + 0.5],
                'r--', linewidth=1, label='Expected')

        ci_95 = np.array([-np.log10(stats.beta.ppf(0.975, np.arange(1, n + 1), n - np.arange(1, n + 1) + 1))])
        ci_05 = np.array([-np.log10(stats.beta.ppf(0.025, np.arange(1, n + 1), n - np.arange(1, n + 1) + 1))])
        ax.fill_between(expected, ci_05.flatten(), ci_95.flatten(),
                         alpha=0.2, color='gray', label='95% CI')

        lambda_gc = np.median(stats.chi2.ppf(1 - df['P'].values, 1)) / stats.chi2.ppf(0.5, 1)
        ax.set_title(f'PWAS QQ Plot ($\\lambda_{{GC}}$ = {lambda_gc:.3f})', fontsize=13)
        ax.set_xlabel('Expected -log10(P)', fontsize=12)
        ax.set_ylabel('Observed -log10(P)', fontsize=12)
        ax.legend(fontsize=10)
        ax.set_aspect('equal')

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] QQ plot 已保存: {output_file}")
        print(f"[INFO] Lambda GC = {lambda_gc:.4f}")

    def plot_volcano(self, result_df, output_file=None, fdr_threshold=0.05):
        if output_file is None:
            output_file = os.path.join(self.output_dir, 'PWAS_Volcano.png')

        df = result_df.dropna(subset=['Z', 'P']).copy()
        df['-log10P'] = -np.log10(df['P'].clip(lower=1e-300))

        fig, ax = plt.subplots(figsize=(10, 8))

        if 'P_FDR' in df.columns:
            sig = df['P_FDR'] < fdr_threshold
            ax.scatter(df.loc[~sig, 'Z'], df.loc[~sig, '-log10P'],
                       s=15, c='gray', alpha=0.5, label='Not significant', edgecolors='none')
            ax.scatter(df.loc[sig, 'Z'], df.loc[sig, '-log10P'],
                       s=30, c='red', alpha=0.8, label=f'FDR < {fdr_threshold}', edgecolors='black', linewidths=0.5)

            for _, row in df[sig].iterrows():
                gene_name = row['ID'].split('.')[-1] if '.' in row['ID'] else row['ID']
                ax.annotate(gene_name, (row['Z'], row['-log10P']),
                            fontsize=7, ha='center', va='bottom',
                            xytext=(0, 5), textcoords='offset points')
        else:
            ax.scatter(df['Z'], df['-log10P'], s=15, c='#1f77b4', alpha=0.6, edgecolors='none')

        ax.axhline(y=-np.log10(0.05 / len(df)), color='red', linestyle='--', linewidth=0.8)
        ax.axvline(x=0, color='gray', linestyle=':', linewidth=0.5)
        ax.set_xlabel('PWAS Z-score', fontsize=12)
        ax.set_ylabel('-log10(P)', fontsize=12)
        ax.set_title('PWAS Volcano Plot', fontsize=14)
        ax.legend(fontsize=10)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] Volcano plot 已保存: {output_file}")

    def compare_weights(self, result_banner, result_rosmap, output_file=None):
        if output_file is None:
            output_file = os.path.join(self.output_dir, 'PWAS_Weight_Comparison.png')

        merged = pd.merge(
            result_banner[['ID', 'Z', 'P']].rename(columns={'Z': 'Z_Banner', 'P': 'P_Banner'}),
            result_rosmap[['ID', 'Z', 'P']].rename(columns={'Z': 'Z_ROSMAP', 'P': 'P_ROSMAP'}),
            on='ID', how='inner'
        )
        merged = merged.dropna()

        if len(merged) == 0:
            print("[WARN] 两个权重集没有共同的基因模型")
            return

        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        axes[0].scatter(merged['Z_Banner'], merged['Z_ROSMAP'],
                        s=10, alpha=0.5, c='#2196F3', edgecolors='none')
        max_z = max(abs(merged['Z_Banner']).max(), abs(merged['Z_ROSMAP']).max())
        axes[0].plot([-max_z, max_z], [-max_z, max_z], 'r--', linewidth=0.8)
        axes[0].set_xlabel('Banner Z-score')
        axes[0].set_ylabel('ROSMAP Z-score')
        axes[0].set_title(f'Z-score Correlation (n={len(merged)})')
        from scipy.stats import pearsonr, spearmanr
        r_p, _ = pearsonr(merged['Z_Banner'], merged['Z_ROSMAP'])
        r_s, _ = spearmanr(merged['Z_Banner'], merged['Z_ROSMAP'])
        axes[0].text(0.05, 0.95, f'Pearson r = {r_p:.3f}\nSpearman ρ = {r_s:.3f}',
                     transform=axes[0].transAxes, fontsize=10, va='top',
                     bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        axes[1].scatter(merged['Z_Banner'], merged['Z_ROSMAP'],
                        s=10, alpha=0.5, c='#2196F3', edgecolors='none')
        axes[1].set_xlabel('Banner Z-score')
        axes[1].set_ylabel('ROSMAP Z-score')
        axes[1].set_title('Z-score Comparison')
        axes[1].plot([-max_z, max_z], [-max_z, max_z], 'r--', linewidth=0.8)

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"[INFO] 权重比较图已保存: {output_file}")

        return merged
