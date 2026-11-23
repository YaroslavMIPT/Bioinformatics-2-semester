import nbformat as nbf
from pathlib import Path

Path("notebooks").mkdir(exist_ok=True)

def md(text):
    return nbf.v4.new_markdown_cell(text.strip())

def code(text):
    return nbf.v4.new_code_cell(text.strip())

# ==============================================================================
# 1. НОУТБУК ПО TCR: nCD4 и nCD8
# ==============================================================================
tcr_cells = [
    md("""
# ДЗ: анализ TCR-репертуаров (nCD4 и nCD8)

Повторяем пайплайн с семинара (где разбирали `nTreg`) для наивных `nCD4` и `nCD8` Т-клеток (цепь TRB). Подбираем нормальные пороги даунсэмплирования под эти выборки, считаем основные метрики, проверяем, как на них влияет глубина, и смотрим, что меняется с возрастом.
"""),
    code("""
import os
from pathlib import Path
import pandas as pd
from pandas.api.types import CategoricalDtype
import matplotlib.pyplot as plt
import seaborn as sns

from repseq import clonosets as cl
from repseq import clone_filter as clf
import repseq
from repseq import stats
stats.AA_PROPS_PATH = os.path.join(os.path.dirname(repseq.__file__), "resourses", "aa_property_table.txt")
from repseq import plot as rsplot
from repseq import beta

PROJECT_ROOT = next(
    p for p in [Path.cwd(), *Path.cwd().parents]
    if (p / "environment.yml").is_file() and (p / "notebooks").is_dir()
)
WORKING_DIR = str(PROJECT_ROOT / "outputs" / "repseq_hw")
os.makedirs(WORKING_DIR, exist_ok=True)

clonosets_ucb_dir = str(PROJECT_ROOT / "data" / "ucb_figshare")
clonosets = cl.find_all_mixcr_clonosets(clonosets_ucb_dir, show_offtarget=False, offtarget_threshold=0.01)

exp_group_order = CategoricalDtype(
    categories=['preterm', 'term', 'late', 'children', 'adults', 'older_adults'],
    ordered=True
)
subset_order = CategoricalDtype(categories=["nCD4", "nTreg", "nCD8"], ordered=True)

metadata = pd.read_csv(os.path.join(clonosets_ucb_dir, "ucb_figshare_TRB_metadata.csv"))[
    ["sample_id", "experimental_group", "donor_id", "subset"]
]
metadata["experimental_group"] = metadata["experimental_group"].astype(exp_group_order)
metadata["subset"] = metadata["subset"].astype(subset_order)
metadata["source"] = metadata["experimental_group"].apply(
    lambda x: "UCB" if x in ['preterm', 'term', 'late'] else "PB"
)

func_filter = clf.Filter(functionality="f", by_umi=True)
print("Всего файлов в папке:", len(clonosets))
"""),
    md("""
## nCD4

Сначала смотрим, сколько образцов в каждой группе и какой у них размер. В `nCD4` всего 30 образцов. Самый маленький — `c2_mnCD4_TRB_1` из группы `children`, у него 5912 функциональных наблюдений (`reads_func`). В группе `children` всего 4 образца, выкидывать её не хочется, поэтому беру `downsample=5500` (с `seed=100`). Потом для проверки устойчивости пересчитаю альфа-разнообразие на `2500`.
"""),
    code("""
clonosets_cd4 = clonosets.merge(metadata).query("subset == 'nCD4'").reset_index(drop=True)
display(clonosets_cd4["experimental_group"].value_counts().sort_index())

stats_cd4 = stats.calc_clonoset_stats(clonosets_cd4, cpu=4)
display(stats_cd4.sort_values("reads_func").head(10)[["sample_id", "clones_func", "reads_func", "clones_nonfunc_freq"]])

rsplot.clonoset_stats(stats_cd4, metadata=metadata, group="experimental_group", properties=["clones_func", "reads_func"], aspect=2)
plt.show()

ds_cd4_main = clf.Filter(functionality="f", downsample=5500, by_umi=True, seed=100)
ds_cd4_alt  = clf.Filter(functionality="f", downsample=2500, by_umi=True, seed=100)
"""),
    md("### Свойства и длины CDR3 (nCD4)"),
    code("""
cdr3_props_cd4 = stats.calc_cdr3_properties(clonosets_cd4, cl_filter=func_filter, cpu=4)
cdr3_props_cd4.to_csv(os.path.join(WORKING_DIR, "cdr3_properties_nCD4.csv"), index=False)
rsplot.cdr3aa_stats(cdr3_props_cd4, metadata=metadata, group="experimental_group", height=3.2, aspect=2)
plt.show()

preterm_late_cd4 = clonosets_cd4.query("experimental_group in ['preterm', 'late']")
cdr3_len_cd4 = stats.cdr3_length_distributions(
    preterm_late_cd4, cl_filter=func_filter, cpu=4,
    count_by_freq=True, seq_type="aa", table="long", zero_fill=True, verbose=False
)
rsplot.cdr3_length_distributions(cdr3_len_cd4, metadata=metadata, group="experimental_group")
plt.show()
"""),
    md("### V-сегменты и конвергенция (nCD4)"),
    code("""
v_usage_cd4 = stats.calc_segment_usage(clonosets_cd4, cl_filter=func_filter, segment="V", cpu=4)
rsplot.segment_usage(v_usage_cd4, metadata=metadata, group="experimental_group")
plt.show()

conv_cd4 = stats.calc_convergence(clonosets_cd4, cl_filter=ds_cd4_main, cpu=4)
rsplot.convergence(conv_cd4, metadata=metadata, group="experimental_group")
plt.show()
"""),
    md("### Альфа-разнообразие, rarefaction и проверка глубины даунсэмпла (nCD4)"),
    code("""
div_cd4_main = stats.calc_diversity_stats(clonosets_cd4, cl_filter=ds_cd4_main, cpu=4)
rsplot.diversity_stats(div_cd4_main, metadata=metadata, group="experimental_group")
plt.show()

rare_cd4 = stats.calc_rarefaction_curve(clonosets_cd4, cl_filter=func_filter, cpu=4)
rsplot.rarefaction_curve(rare_cd4)
plt.show()

# смотрим, что будет при downsample=2500
div_cd4_alt = stats.calc_diversity_stats(clonosets_cd4, cl_filter=ds_cd4_alt, cpu=4)
rsplot.diversity_stats(div_cd4_alt, metadata=metadata, group="experimental_group")
plt.show()
"""),
    md("### Бета-разнообразие (F2) и MDS (nCD4)"),
    code("""
beta_cd4 = beta.metrics(clonosets_cd4, cl_filter=ds_cd4_main, cpu=4)
rsplot.beta_table(beta_cd4)
plt.show()

mds_cd4 = beta.mds(beta_cd4, metric="f2")
rsplot.beta_mds(mds_cd4, metadata=metadata, group="experimental_group")
plt.show()
"""),
    md("""
## nCD8

Теперь то же самое для `nCD8`. Тут 27 образцов, самые маленькие попали в `older_adults`, минимум у `o7_nCD8_TRB_1` — 6368 `reads_func`. Чтобы не терять образцы из пожилой группы, ставим `downsample=6000` (`seed=100`), а для проверки чувствительности сравним с `3000`.
"""),
    code("""
clonosets_cd8 = clonosets.merge(metadata).query("subset == 'nCD8'").reset_index(drop=True)
display(clonosets_cd8["experimental_group"].value_counts().sort_index())

stats_cd8 = stats.calc_clonoset_stats(clonosets_cd8, cpu=4)
display(stats_cd8.sort_values("reads_func").head(10)[["sample_id", "clones_func", "reads_func", "clones_nonfunc_freq"]])

rsplot.clonoset_stats(stats_cd8, metadata=metadata, group="experimental_group", properties=["clones_func", "reads_func"], aspect=2)
plt.show()

ds_cd8_main = clf.Filter(functionality="f", downsample=6000, by_umi=True, seed=100)
ds_cd8_alt  = clf.Filter(functionality="f", downsample=3000, by_umi=True, seed=100)
"""),
    md("### Свойства и длины CDR3 (nCD8)"),
    code("""
cdr3_props_cd8 = stats.calc_cdr3_properties(clonosets_cd8, cl_filter=func_filter, cpu=4)
cdr3_props_cd8.to_csv(os.path.join(WORKING_DIR, "cdr3_properties_nCD8.csv"), index=False)
rsplot.cdr3aa_stats(cdr3_props_cd8, metadata=metadata, group="experimental_group", height=3.2, aspect=2)
plt.show()

preterm_late_cd8 = clonosets_cd8.query("experimental_group in ['preterm', 'late']")
cdr3_len_cd8 = stats.cdr3_length_distributions(
    preterm_late_cd8, cl_filter=func_filter, cpu=4,
    count_by_freq=True, seq_type="aa", table="long", zero_fill=True, verbose=False
)
rsplot.cdr3_length_distributions(cdr3_len_cd8, metadata=metadata, group="experimental_group")
plt.show()
"""),
    md("### V-сегменты и конвергенция (nCD8)"),
    code("""
v_usage_cd8 = stats.calc_segment_usage(clonosets_cd8, cl_filter=func_filter, segment="V", cpu=4)
rsplot.segment_usage(v_usage_cd8, metadata=metadata, group="experimental_group")
plt.show()

conv_cd8 = stats.calc_convergence(clonosets_cd8, cl_filter=ds_cd8_main, cpu=4)
rsplot.convergence(conv_cd8, metadata=metadata, group="experimental_group")
plt.show()
"""),
    md("### Альфа-разнообразие, rarefaction и проверка глубины даунсэмпла (nCD8)"),
    code("""
div_cd8_main = stats.calc_diversity_stats(clonosets_cd8, cl_filter=ds_cd8_main, cpu=4)
rsplot.diversity_stats(div_cd8_main, metadata=metadata, group="experimental_group")
plt.show()

rare_cd8 = stats.calc_rarefaction_curve(clonosets_cd8, cl_filter=func_filter, cpu=4)
rsplot.rarefaction_curve(rare_cd8)
plt.show()

# смотрим, что будет при downsample=3000
div_cd8_alt = stats.calc_diversity_stats(clonosets_cd8, cl_filter=ds_cd8_alt, cpu=4)
rsplot.diversity_stats(div_cd8_alt, metadata=metadata, group="experimental_group")
plt.show()
"""),
    md("### Бета-разнообразие (F2) и MDS (nCD8)"),
    code("""
beta_cd8 = beta.metrics(clonosets_cd8, cl_filter=ds_cd8_main, cpu=4)
rsplot.beta_table(beta_cd8)
plt.show()

mds_cd8 = beta.mds(beta_cd8, metric="f2")
rsplot.beta_mds(mds_cd8, metadata=metadata, group="experimental_group")
plt.show()
"""),
    md("""
## Выводы по nCD4 и nCD8

В целом картина по возрасту в `nCD4` и `nCD8` похожа на то, что было в `nTreg`, но со своими нюансами.

У недоношенных (`preterm`) в обеих популяциях CDR3 заметно короче, среднее число вставок меньше, а доля последовательностей вообще без вставок (`zero_insert_freq`) выше, чем у доношенных и тем более у детей и взрослых. Это хорошо ложится на биологию: на ранних сроках гестации активность TdT ещё низкая, поэтому случайных N-нуклеотидов вставляется мало.

Из-за этого же в пуповинной крови (особенно в `preterm`) выше конвергенция и попарное сходство F2 между разными донорами — без длинных случайных вставок чаще независимо собираются одинаковые последовательности. С возрастом репертуары становятся более «индивидуальными», и на MDS группы пуповинной и периферической крови расходятся.

Главное отличие между `nCD4` и `nCD8` видно по альфа-разнообразию у пожилых (`older_adults`). Если наивные `nCD4` держат довольно высокое разнообразие почти всю жизнь, то в `nCD8` у `older_adults` разнообразие проседает гораздо сильнее, а клональность растёт — сказывается инволюция тимуса и гомеостатическая экспансия отдельных CD8-клонов.

Когда уменьшили `downsample` примерно вдвое (5500 → 2500 для `nCD4` и 6000 → 3000 для `nCD8`), абсолютные числа наблюдаемых клонотипов и Chao1, естественно, упали, потому что в меньшую выборку попадает меньше редких клонов. Но сами соотношения между возрастными группами не сломались — тренды остались теми же.
""")
]

nb_tcr = nbf.v4.new_notebook(cells=tcr_cells)
nbf.write(nb_tcr, "notebooks/UCB_homework_TCR_CD4_CD8.ipynb")

# ==============================================================================
# 2. НОУТБУК ПО BCR (MiXCR -> repseq)
# ==============================================================================
bcr_cells = [
    md("""
# ДЗ по BCR: обработка в MiXCR и анализ в repseq

Исходные парные риды `DChu_T7_pbmc_d7COVID2_2` (PBMC на 7-й день после COVID-19) прогнал через MiXCR 4.7.0 с пресетом `generic-amplicon` (`--species hsa --rna --floating-left-alignment-boundary --floating-right-alignment-boundary C --assemble-clonotypes-by CDR3`).
По QC всё прошло чисто: выровнялось 93.08% ридов, off-target всего 0.1%, собралась тяжёлая цепь `IGH` (49 079 клонотипов, из них 46 795 функциональных, 760 664 функциональных рида).

Так как ноутбук с семинара был написан под 30 образцов TCR, напрямую его не перенесёшь: здесь всего один образец, поэтому межгрупповые сравнения, F2 и MDS смысла не имеют. Ниже оставил те метрики из `repseq`, которые работают на одиночном репертуаре, и отдельно посмотрел специфичные для BCR вещи — топ-клоны, V-сегменты и изотипы.
"""),
    code("""
import os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from repseq import clonosets as cl
from repseq import clone_filter as clf
import repseq
from repseq import stats
stats.AA_PROPS_PATH = os.path.join(os.path.dirname(repseq.__file__), "resourses", "aa_property_table.txt")
from repseq import plot as rsplot

PROJECT_ROOT = next(
    p for p in [Path.cwd(), *Path.cwd().parents]
    if (p / "environment.yml").is_file() and (p / "notebooks").is_dir()
)
bcr_dir = str(PROJECT_ROOT / "outputs" / "bcr_mixcr")
bcr_clonosets = cl.find_all_mixcr_clonosets(bcr_dir, show_offtarget=False, offtarget_threshold=0.01)
display(bcr_clonosets)

# UMI-колонки тут нет, считаем по ридам
func_filter = clf.Filter(functionality="f", by_umi=False)
ds_filter   = clf.Filter(functionality="f", downsample=50000, by_umi=False, seed=100)
"""),
    md("## Статистика клоносета"),
    code("""
bcr_stats = stats.calc_clonoset_stats(bcr_clonosets, cpu=2)
display(bcr_stats[["sample_id", "chain", "clones", "clones_func", "clones_nonfunc_freq", "reads", "reads_func"]])
"""),
    md("## Свойства и распределение длин CDR3"),
    code("""
bcr_cdr3_props = stats.calc_cdr3_properties(bcr_clonosets, cl_filter=func_filter, cpu=2)
display(bcr_cdr3_props[["sample_id", "mean_cdr3nt_len", "mean_insert_size", "zero_insert_freq", "cdr3_full_hydropathy", "cdr3_full_charge"]])

bcr_len_dist = stats.cdr3_length_distributions(
    bcr_clonosets, cl_filter=func_filter, cpu=2,
    count_by_freq=True, seq_type="aa", table="long", zero_fill=True, verbose=False
)
rsplot.cdr3_length_distributions(bcr_len_dist)
plt.title("Длины CDR3 (aa) — BCR IGH")
plt.show()
"""),
    md("## Разнообразие, конвергенция и rarefaction"),
    code("""
bcr_div = stats.calc_diversity_stats(bcr_clonosets, cl_filter=ds_filter, cpu=2)
display(bcr_div)

bcr_conv = stats.calc_convergence(bcr_clonosets, cl_filter=ds_filter, cpu=2)
display(bcr_conv)

bcr_rare = stats.calc_rarefaction_curve(bcr_clonosets, cl_filter=func_filter, cpu=2)
rsplot.rarefaction_curve(bcr_rare)
plt.title("Rarefaction — BCR IGH")
plt.show()
"""),
    md("## Топ-клонотипы, V-гены и изотипы"),
    code("""
raw_bcr = pd.read_csv(bcr_clonosets.loc[0, "filename"], sep="\\t")
raw_bcr["topV"] = raw_bcr["allVHitsWithScore"].astype(str).str.split("*").str[0]

print("Топ-10 клонотипов IGH:")
cols = [c for c in ["cloneId", "readCount", "readFraction", "aaSeqCDR3", "allVHitsWithScore", "allCHitsWithScore"] if c in raw_bcr.columns]
display(raw_bcr[cols].head(10))

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

top_v = raw_bcr.groupby("topV")["readFraction"].sum().sort_values(ascending=False).head(15)
top_v.plot(kind="bar", ax=axes[0], color="steelblue")
axes[0].set_title("Топ-15 V-сегментов (по доле ридов)")
axes[0].set_ylabel("readFraction")

if "allCHitsWithScore" in raw_bcr.columns:
    raw_bcr["isotype"] = raw_bcr["allCHitsWithScore"].astype(str).str.split("*").str[0].str.replace("nan", "Unknown")
    iso = raw_bcr.groupby("isotype")["readFraction"].sum().sort_values(ascending=False)
    iso.plot(kind="bar", ax=axes[1], color="darkorange")
    axes[1].set_title("Изотипы (C-сегменты)")
    axes[1].set_ylabel("readFraction")
else:
    axes[1].text(0.5, 0.5, "C-сегменты не экспортированы", ha="center")

plt.tight_layout()
plt.show()
"""),
    md("""
## Что видно по BCR

Образец просеквенирован глубоко (760 тыс. функциональных ридов, почти 47 тыс. клонотипов `IGH`), при этом доля нефункциональных перестроек небольшая.

Первое, что бросается в глаза при сравнении с TCR из предыдущего ноутбука — распределение длин CDR3 у `IGH` заметно шире и сдвинуто в сторону более длинных последовательностей, что нормально для тяжёлых цепей антител. Второе — сильная клональная экспансия: один только самый крупный клон забирает больше 1.6% всех ридов, что хорошо согласуется с 7-м днём после COVID-19 (пик ответа плазмабластов).

Ну и главное ограничение такого переноса TCR-пайплайна на BCR: сборка и подсчёт по точному совпадению CDR3 не учитывают соматическую гипермутацию (SHM). Потомки одного и того же B-клеточного клона с точечными мутациями здесь считаются как разные клонотипы, что искусственно раздувает разнообразие и не даёт увидеть структуру клональных деревьев и переключение классов. Для полноценного анализа BCR нужно собирать полные VDJ-последовательности и строить деревья гипермутаций.
""")
]

nb_bcr = nbf.v4.new_notebook(cells=bcr_cells)
nbf.write(nb_bcr, "notebooks/BCR_homework.ipynb")
print("Обновлённые ноутбуки созданы.")
