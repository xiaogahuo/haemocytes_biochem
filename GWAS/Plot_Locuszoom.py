
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from SUtils import bed_path, exp_name

def plot_locuszoom(df, file_name, chr_num, start_pos, end_pos, sig_level=5e-8, window_size=None):
    """
    生成指定区域的 LocusZoom 图
    :param file_name: 数据文件名, 如 '0_VS_1'
    :param chr_num: 染色体编号 (1-22)
    :param start_pos: 区域起始位置 (bp)
    :param end_pos: 区域结束位置 (bp)
    :param sig_level: 显著性阈值线, 默认 5e-8
    :param window_size: 如果指定, 则以 lead SNP 为中心扩展窗口 (bp)
    """
    print("columns: {}".format(df.columns.tolist()))
    print("total snps: {}".format(df.shape[0]))

    df_region = df[(df['CHR'] == chr_num) & (df['BP'] >= start_pos) & (df['BP'] <= end_pos)].copy()
    print("region snps: {}".format(df_region.shape[0]))

    if df_region.empty:
        print("WARNING: No SNPs found in chr{}:{}-{}".format(chr_num, start_pos, end_pos))
        return

    df_region['-log10P'] = -np.log10(df_region['P'].clip(lower=1e-300))

    lead_idx = df_region['P'].idxmin()
    lead_snp = df_region.loc[lead_idx, 'SNP']
    lead_bp = df_region.loc[lead_idx, 'BP']
    lead_p = df_region.loc[lead_idx, 'P']
    print("lead SNP: {}, pos: {}, P: {:.2e}".format(lead_snp, lead_bp, lead_p))

    if window_size is not None:
        half_w = window_size // 2
        start_pos = max(1, lead_bp - half_w)
        end_pos = lead_bp + half_w
        df_region = df[(df['CHR'] == chr_num) & (df['BP'] >= start_pos) & (df['BP'] <= end_pos)].copy()
        df_region['-log10P'] = -np.log10(df_region['P'].clip(lower=1e-300))
        print("after window adjustment, region snps: {}".format(df_region.shape[0]))

    try:
        recomb_df = _load_recomb_data(chr_num, start_pos, end_pos)
    except Exception as e:
        print("recombination data not available: {}".format(e))
        recomb_df = None

    fig, axes = plt.subplots(2 if recomb_df is not None else 1, 1,
                             figsize=(12, 6 if recomb_df is None else 8),
                             gridspec_kw={'height_ratios': [4, 1]} if recomb_df is not None else None,
                             dpi=300)
    if recomb_df is None:
        axes = [axes]

    ax_manh = axes[0]
    scatter = ax_manh.scatter(df_region['BP'], df_region['-log10P'],
                              c=df_region['P'], cmap='RdYlGn_r', vmin=0, vmax=1,
                              s=30, edgecolors='black', linewidths=0.3, alpha=0.8, zorder=3)

    ax_manh.axhline(y=-np.log10(sig_level), color='red', linestyle='--', linewidth=1, alpha=0.7, zorder=1)
    ax_manh.axhline(y=-np.log10(1e-5), color='blue', linestyle=':', linewidth=0.8, alpha=0.5, zorder=1)

    if lead_p < sig_level:
        ax_manh.annotate(lead_snp, xy=(lead_bp, -np.log10(lead_p)),
                         xytext=(lead_bp + (end_pos - start_pos) * 0.05, -np.log10(lead_p) + 0.5),
                         fontsize=9, fontweight='bold', color='darkred',
                         arrowprops=dict(arrowstyle='->', color='darkred', lw=1.2), zorder=5)

    ax_manh.set_xlabel('Chromosome {} Position (bp)'.format(chr_num), fontsize=11)
    ax_manh.set_ylabel('$-log_{10}(P)$', fontsize=11)
    ax_manh.set_title('LocusZoom: {} (chr{}:{}-{})'.format(file_name, chr_num, start_pos, end_pos), fontsize=12)
    ax_manh.set_xlim(start_pos, end_pos)
    ax_manh.grid(True, linestyle='--', alpha=0.3, zorder=0)

    cbar = plt.colorbar(scatter, ax=ax_manh, shrink=0.6, pad=0.02)
    cbar.set_label('$r^2$ with {}'.format(lead_snp), fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    if recomb_df is not None and len(axes) > 1:
        ax_rec = axes[1]
        ax_rec.fill_between(recomb_df['BP'], recomb_df['Rate'], color='#888888', alpha=0.6)
        ax_rec.set_xlabel('Chromosome {} Position (bp)'.format(chr_num), fontsize=11)
        ax_rec.set_ylabel('Recombination\nRate (cM/Mb)', fontsize=9)
        ax_rec.set_xlim(start_pos, end_pos)
        ax_rec.grid(True, linestyle='--', alpha=0.3)

    plt.tight_layout()
    out_dir = bed_path + "{}/{}/locuszoom/".format(exp_name, file_name)
    os.makedirs(out_dir, exist_ok=True)
    out_path = out_dir + "LocusZoom_{}_chr{}_{}_{}.png".format(file_name, chr_num, start_pos, end_pos)
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.show()
    plt.clear()
    print("LocusZoom plot saved to: {}".format(out_path))


def plot_locuszoom_from_lead(file_name, lead_snp, window_size=500000, sig_level=5e-8):
    """
    以指定的 lead SNP 为中心, 自动生成 LocusZoom 图
    :param file_name: 数据文件名
    :param lead_snp: lead SNP 名称
    :param window_size: 窗口大小 (bp), 默认 500kb
    :param sig_level: 显著性阈值
    """
    summary_path = bed_path + "{}/{}/{}_gwas_summary.txt.gz".format(exp_name, file_name, file_name)
    df = pd.read_csv(summary_path, sep=r'\s+', engine='python')

    snp_row = df[df['SNP'] == lead_snp]
    if snp_row.empty:
        print("ERROR: SNP {} not found in data".format(lead_snp))
        return

    chr_num = int(snp_row.iloc[0]['CHR'])
    bp_pos = int(snp_row.iloc[0]['BP'])
    half_w = window_size // 2
    start_pos = max(1, bp_pos - half_w)
    end_pos = bp_pos + half_w

    print("plotting LocusZoom for {} at chr{}:{} (window: {}kb)".format(lead_snp, chr_num, bp_pos, window_size // 1000))
    plot_locuszoom(df, file_name, chr_num, start_pos, end_pos, sig_level=sig_level)


def plot_locuszoom_batch(file_name, sig_level=5e-8, window_size=500000, top_n=10):
    """
    批量生成 top N 独立位点的 LocusZoom 图
    :param file_name: 数据文件名
    :param sig_level: 显著性阈值
    :param window_size: 窗口大小 (bp)
    :param top_n: 要绘制的独立位点数量
    """
    summary_path = bed_path + "{}/{}/{}_gwas_summary.txt.gz".format(exp_name, file_name, file_name)
    df = pd.read_csv(summary_path, sep=r'\s+', engine='python')
    df_sig = df[df['P'] < sig_level].sort_values('P')

    if df_sig.empty:
        print("No significant SNPs found at threshold {}".format(sig_level))
        return

    lead_snps = []
    plotted_regions = []
    for _, row in df_sig.iterrows():
        chr_num, bp = int(row['CHR']), int(row['BP'])
        is_independent = True
        for p_chr, p_start, p_end in plotted_regions:
            if chr_num == p_chr and p_start <= bp <= p_end:
                is_independent = False
                break
        if is_independent:
            half_w = window_size // 2
            plotted_regions.append((chr_num, max(1, bp - half_w), bp + half_w))
            lead_snps.append(row['SNP'])
            if len(lead_snps) >= top_n:
                break

    print("found {} independent loci to plot".format(len(lead_snps)))
    for i, snp in enumerate(lead_snps):
        print("[{}/{}] plotting {}...".format(i + 1, len(lead_snps), snp))
        plot_locuszoom_from_lead(file_name, snp, window_size=window_size, sig_level=sig_level)


def _load_recomb_data(chr_num, start_pos, end_pos):
    """加载重组率数据 (如果存在)"""
    recomb_path = bed_path + "recomb/recomb_rate_chr{}.txt".format(chr_num)
    if not os.path.exists(recomb_path):
        return None
    df = pd.read_csv(recomb_path, sep=r'\s+', engine='python')
    df_region = df[(df['BP'] >= start_pos) & (df['BP'] <= end_pos)]
    return df_region


if __name__ == '__main__':
    # === 方式1: 指定染色体区域绘制 ===
    # plot_locuszoom(f_name, chr_num=6, start_pos=25000000, end_pos=35000000)

    # === 方式2: 以指定 lead SNP 为中心绘制 ===
    # plot_locuszoom_from_lead("0_VS_1", lead_snp="rs199504", window_size=500000)
    # plot_locuszoom_from_lead("0_VS_2", lead_snp="rs6936540", window_size=500000)

    # === 方式3: 批量绘制 top N 独立位点 ===
    plot_locuszoom_batch("0_VS_12", sig_level=5e-8, window_size=500000, top_n=2)