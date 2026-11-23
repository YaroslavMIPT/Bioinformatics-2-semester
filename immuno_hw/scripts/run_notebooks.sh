#!/bin/bash
#SBATCH --job-name=run_hw_nbs
#SBATCH --partition=AMD9554-common
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=outputs/run_hw_nbs_%j.log

source /home/STUDY/FBMF/studfbmf02_13/miniforge3/etc/profile.d/conda.sh
conda activate immune-mixcr

python scripts/build_notebooks.py

echo ">>> Executing TCR notebook..."
jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.timeout=1200 notebooks/UCB_homework_TCR_CD4_CD8.ipynb

echo ">>> Executing BCR notebook..."
jupyter nbconvert --to notebook --execute --inplace --allow-errors --ExecutePreprocessor.timeout=1200 notebooks/BCR_homework.ipynb

echo ">>> Checking for any cell errors..."
python -c "
import nbformat
for path in ['notebooks/UCB_homework_TCR_CD4_CD8.ipynb', 'notebooks/BCR_homework.ipynb']:
    nb = nbformat.read(path, as_version=4)
    errs = []
    for i, c in enumerate(nb.cells):
        for out in c.get('outputs', []):
            if out.get('output_type') == 'error':
                errs.append((i, out.get('ename'), out.get('evalue')))
    print(path, '-> errors:', len(errs))
    for e in errs:
        print('  cell', e)
"
