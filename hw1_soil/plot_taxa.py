import os
import pandas as pd
import matplotlib.pyplot as plt
import qiime2

BASE_DIR = '/home/STUDY/FBMF/studfbmf02_13/metagenome_soil'
FIG_DIR = os.path.join(BASE_DIR, 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

# 1. Загрузка данных
meta = pd.read_csv('/home/STUDY/FBMF/bioinformatics/metagenomes/soil_hw/soil_metadata_full.tsv', 
                   sep='\t', index_col=0, skiprows=[1])
table = qiime2.Artifact.load(os.path.join(BASE_DIR, 'qza/soil_ASV_table.qza')).view(pd.DataFrame)
taxa = qiime2.Artifact.load(os.path.join(BASE_DIR, 'qza/taxonomy.qza')).view(pd.Series)

# 2. Агрегация по Phylum
def extract_phylum(taxon_str):
    for part in taxon_str.split(';'):
        part = part.strip()
        if part.startswith('p__'):
            name = part.replace('p__', '')
            return name if name else 'Unassigned'
    return 'Unassigned'

phyla = taxa.apply(extract_phylum)
phylum_table = table.T.groupby(phyla).sum().T
phylum_rel = phylum_table.div(phylum_table.sum(axis=1), axis=0) * 100

# Выделение топ-7 филумов, остальные в 'Other'
top_phyla = phylum_rel.mean().sort_values(ascending=False).head(7).index.tolist()
plot_df = phylum_rel[top_phyla].copy()
plot_df['Other'] = phylum_rel.drop(columns=top_phyla).sum(axis=1)

# Объединение с метаданными и усреднение по contamination_level
plot_df['contamination_level'] = meta['contamination_level']
grouped = plot_df.groupby('contamination_level').mean().reindex(['k0', 'k5', 'k25'])

# 3. Отрисовка stacked barplot
plt.figure(figsize=(8, 6))
colors = ['#1f77b4', '#aec7e8', '#ff7f0e', '#2ca02c', '#98df8a', '#d62728', '#9467bd', '#c7c7c7']
ax = grouped.plot(kind='bar', stacked=True, figsize=(8, 6), color=colors, edgecolor='black', width=0.6)

plt.title('Таксономический состав сообщества (Phylum)', weight='bold', pad=12)
plt.xlabel('Уровень загрязнения', weight='bold')
plt.ylabel('Относительное обилие (%)', weight='bold')
plt.xticks(ticks=[0, 1, 2], labels=['Контроль (k0)', 'Низкое (k5)', 'Высокое (k25)'], rotation=0)
plt.ylim(0, 100)
plt.legend(bbox_to_anchor=(1.04, 1), loc='upper left', title='Филум', frameon=True)
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()

out_png = os.path.join(FIG_DIR, 'taxa_barplot.png')
out_pdf = os.path.join(FIG_DIR, 'taxa_barplot.pdf')
plt.savefig(out_png, dpi=300)
plt.savefig(out_pdf)
plt.close()
print(f'Сохранен: {out_png}')
