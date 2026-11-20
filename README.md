# Which nucleus features separate malignant from benign breast masses?

Capstone analysis of a real public dataset, with multiple-testing control and
a cautious written interpretation.

## Biological question

Do measurements of cell-nucleus shape and size, taken from digitised images of
fine-needle aspirates of breast masses, differ between malignant and benign
masses, and which of the 30 measured features differ most?

## Dataset

Wisconsin Diagnostic Breast Cancer (WDBC), UCI Machine Learning Repository
(Wolberg, Street and Mangasarian). 569 samples: 212 malignant, 357 benign.
Each sample has 30 numeric features (10 nucleus measurements such as radius,
texture and concavity, each summarised as mean, standard error and "worst").
I use the copy bundled with scikit-learn, so no download is needed. The script
saves it as `data/wdbc.csv`.

## Files

| File | What it is |
|---|---|
| `analysis.py` | The whole analysis |
| `requirements.txt` | Python packages needed |
| `.github/workflows/analysis.yml` | Runs the analysis on GitHub (manual button) |
| `data/wdbc.csv` | The dataset used (written by the script) |
| `results/feature_tests.tsv` | Per-feature statistics, p-values and corrected values |
| `results/summary.json` | Headline numbers, including the permutation control |
| `results/effect_vs_significance.png` | Effect size against significance |
| `results/feature_correlations.png` | Feature correlation heatmap |
| `INTERPRETATION.md` | Written interpretation, limitations and confounders |

## How to run

On GitHub: **Actions** tab -> **Run analysis** -> **Run workflow**. Results are
committed to `data/` and `results/`.

On your computer:

    pip install -r requirements.txt
    python3 analysis.py

The random seed is fixed (`SEED = 1`) so the permutation control is repeatable.

## Methods in brief

1. For each of the 30 features, a two-sided Mann-Whitney U test compares
   malignant with benign samples.
2. Effect size is the rank-biserial correlation (+1 means malignant values are
   always higher; 0 means no difference).
3. Multiple testing: p-values are adjusted with Benjamini-Hochberg (FDR 5%,
   main method) and Bonferroni (stricter, for comparison).
4. Permutation control: diagnosis labels are shuffled 200 times and the same
   pipeline is rerun, to show how many false positives appear with and without
   correction when there is no real signal.
5. Feature-feature correlations are checked, because correlated features mean
   the 30 tests are not independent.

## Decisions

- Mann-Whitney U instead of a t-test, because many features are skewed.
- BH as the main correction because the tests are many and correlated, and FDR
  control is the usual choice for feature screening. Bonferroni is shown as a
  stricter check.
- I report effect sizes next to p-values, because with 569 samples even tiny
  differences become significant.
- No classifier was built: the question is about which features differ, not
  about diagnosis.

See `INTERPRETATION.md` for results, limitations and confounders.
