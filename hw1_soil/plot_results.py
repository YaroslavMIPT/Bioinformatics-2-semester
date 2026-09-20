import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skbio import OrdinationResults
import qiime2

BASE_DIR = '/home/STUDY/FBMF/studfbmf02_13/metagenome_soil'
FIG_DIR = os.path.join(BASE_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

sns.set_theme(style='ticks', font_scale=1.1)

# 1. Метаданные
meta_file = '/home/STUDY/FBMF/bioinformatics/metagenomes/soil_hw/soil_metadata_full.tsv'
meta_df = pd.read_csv(meta_file, sep='\t', index_col=0, skiprows=[1])

contam_order = ['k0', 'k5', 'k25']
day_order = ['d3', 'd90', 'd180', 'd360']
palette = {'k0': '#2ca02c', 'k5': '#ff7f0e', 'k25': '#d62728'}
marker_dict = {'d3': 'o', 'd90': 's', 'd180': '^', 'd360': 'D'}

# 2. Построение PCoA напрямую из .qza
def plot_pcoa_from_qza(qza_path, title, filename):
    pcoa_art = qiime2.Artifact.load(qza_path)
    ord_res = pcoa_art.view(OrdinationResults)
    
    coords = ord_res.samples.iloc[:, :2].copy()
    coords.columns = ['PC1', 'PC2']
    var1 = ord_res.proportion_explained[0] * 100
    var2 = ord_res.proportion_explained[1] * 100

    merged = coords.join(meta_df, how='inner')

    plt.figure(figsize=(7.5, 6))
    sns.scatterplot(
        data=merged,
        x='PC1', y='PC2',
        hue='contamination_level',
        style='day_after_contamination',
        hue_order=contam_order,
        style_order=day_order,
        palette=palette,
        markers=marker_dict,
        s=130,
        alpha=0.9,
        edgecolor='k',
        linewidth=0.8
    )
    plt.title(title, weight='bold', pad=12)
    plt.xlabel(f'PC1 ({var1:.1f}%)', weight='bold')
    plt.ylabel(f'PC2 ({var2:.1f}%)', weight='bold')
    plt.grid(True, linestyle='--', alpha=0.4)
    plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left', frameon=True)
    plt.tight_layout()
    
    out_png = os.path.join(FIG_DIR, f'{filename}.png')
    out_pdf = os.path.join(FIG_DIR, f'{filename}.pdf')
    plt.savefig(out_png, dpi=300)
    plt.savefig(out_pdf)
    plt.close()
    print(f'Сохранен: {out_png}')

plot_pcoa_from_qza(
    os.path.join(BASE_DIR, 'core-metrics-results/bray_curtis_pcoa_results.qza'),
    'PCoA — Bray-Curtis (PERMANOVA p = 0.001)',
    'pcoa_bray_curtis'
)

plot_pcoa_from_qza(
    os.path.join(BASE_DIR, 'core-metrics-results/unweighted_unifrac_pcoa_results.qza'),
    'PCoA — Unweighted UniFrac (PERMANOVA p = 0.009)',
    'pcoa_unweighted_unifrac'
)

# 3. Боксплоты альфа-разнообразия
shannon = qiime2.Artifact.load(os.path.join(BASE_DIR, 'core-metrics-results/shannon_vector.qza')).view(pd.Series)
faith = qiime2.Artifact.load(os.path.join(BASE_DIR, 'core-metrics-results/faith_pd_vector.qza')).view(pd.Series)

alpha_df = meta_df.copy()
alpha_df['Shannon'] = shannon
alpha_df['Faith_PD'] = faith

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Shannon по дозам
sns.boxplot(data=alpha_df, x='contamination_level', y='Shannon', order=contam_order, 
            palette=palette, ax=axes[0, 0], boxprops=dict(alpha=0.8))
sns.stripplot(data=alpha_df, x='contamination_level', y='Shannon', order=contam_order, 
              color='black', size=6, jitter=0.2, ax=axes[0, 0])
axes[0, 0].set_title('Shannon: уровень загрязнения (p = 0.008)', weight='bold')
axes[0, 0].set_xlabel('Уровень загрязнения')
axes[0, 0].set_ylabel('Индекс Шеннона')

# Shannon по дням
sns.boxplot(data=alpha_df, x='day_after_contamination', y='Shannon', order=day_order, 
            color='#aec7e8', ax=axes[0, 1], boxprops=dict(alpha=0.8))
sns.stripplot(data=alpha_df, x='day_after_contamination', y='Shannon', order=day_order, 
              color='black', size=6, jitter=0.2, ax=axes[0, 1])
axes[0, 1].set_title('Shannon: динамика во времени (p = 0.0003)', weight='bold')
axes[0, 1].set_xlabel('Сутки после внесения')
axes[0, 1].set_ylabel('Индекс Шеннона')

# Faith PD по дозам
sns.boxplot(data=alpha_df, x='contamination_level', y='Faith_PD', order=contam_order, 
            palette=palette, ax=axes[1, 0], boxprops=dict(alpha=0.8))
sns.stripplot(data=alpha_df, x='contamination_level', y='Faith_PD', order=contam_order, 
              color='black', size=6, jitter=0.2, ax=axes[1, 0])
axes[1, 0].set_title("Faith's PD: уровень загрязнения (p = 0.057)", weight='bold')
axes[1, 0].set_xlabel('Уровень загрязнения')
axes[1, 0].set_ylabel("Faith's PD")

# Faith PD по дням
sns.boxplot(data=alpha_df, x='day_after_contamination', y='Faith_PD', order=day_order, 
            color='#98df8a', ax=axes[1, 1], boxprops=dict(alpha=0.8))
sns.stripplot(data=alpha_df, x='day_after_contamination', y='Faith_PD', order=day_order, 
              color='black', size=6, jitter=0.2, ax=axes[1, 1])
axes[1, 1].set_title("Faith's PD: динамика во времени (p = 0.0003)", weight='bold')
axes[1, 1].set_xlabel('Сутки после внесения')
axes[1, 1].set_ylabel("Faith's PD")

plt.tight_layout()
out_box_png = os.path.join(FIG_DIR, 'alpha_diversity_boxplots.png')
out_box_pdf = os.path.join(FIG_DIR, 'alpha_diversity_boxplots.pdf')
plt.savefig(out_box_png, dpi=300)
plt.savefig(out_box_pdf)
plt.close()
print(f'Сохранен: {out_box_png}')
