# Interpretation

## Question

Do nucleus size and shape measurements differ between malignant and benign
breast masses in the WDBC dataset (212 malignant, 357 benign), and which
features differ most?

## What I found

- **27 of 30 features differ** between malignant and benign samples at a 5%
  false discovery rate (Benjamini-Hochberg). Bonferroni keeps 26.
- **Differences are largest for size and concavity features.** The top five are
  worst perimeter, worst radius, worst area, worst concave points and mean
  concave points, with rank-biserial effect sizes of 0.93 to 0.95. An effect
  size of 0.95 corresponds to roughly a 97.5% chance that a random malignant
  sample has a higher value than a random benign one. Malignant samples had
  larger and more irregular nuclei.
- **Three features show no clear difference:** smoothness error, mean fractal
  dimension and texture error (q = 0.23 to 0.64).
- **One result is borderline:** symmetry error passes BH (q = 0.031) but fails
  Bonferroni, and its effect is small (-0.11), in the opposite direction
  (higher in benign). I treat it as unconfirmed.

## What the multiple-testing control did

On the real data, correction changed little (27 significant before and after
BH), because the strongest signals are far below any threshold. The correction
matters when there is no real signal. In the permutation control, with
diagnosis labels shuffled 200 times:

- uncorrected testing gave on average 1.6 false positives per run, and 94 of 200 runs had at
  least one;
- BH gave on average 0.14 false positives per run, and 6 of 200 runs (3%) had at least one.

So uncorrected testing at 5% reports "findings" in about half of pure-noise
datasets, while BH keeps that near or below 5%. BH controls the expected share
of false discoveries, not the chance of any: with many real signals, some of
the 27 could still be false.

## Limitations and confounders

- **The 30 tests are not independent.** 26 of 435 feature pairs have a
  Spearman correlation of 0.9 or more. Radius, perimeter and area are
  geometrically linked, and the mean, error and "worst" versions of a feature
  come from the same nuclei. The "27 significant features" therefore reflect a
  much smaller number of independent signals, mainly overall nucleus size and
  contour irregularity. BH is valid under positive dependence, which is
  plausible here, but I did not prove it for this data.
- **Possible repeated samples per patient.** The dataset has no patient
  identifiers, so I cannot check whether one patient contributed several
  samples. If so, samples are not independent and the p-values are too
  optimistic.
- **Circularity risk.** Nuclear size and shape are what cytologists look at
  when judging malignancy. If the diagnosis labels relied partly on such
  appearance, the differences are partly built in. I did not verify how
  diagnoses were confirmed, so I cannot rule this out.
- **Unmeasured confounders.** Tumour size and grade, patient age, menopausal
  status, sampling site and sample preparation or imaging batch are not in the
  data. Any of these could differ between benign and malignant groups and also
  affect nucleus measurements, so the differences are associations, not shown
  to be caused by malignancy alone.
- **Selection and generalisability.** The samples come from one institution and
  were referred for aspiration, so 37% malignant is not a population rate.
  Results may not hold for other populations, scanners or staining methods.
- **Statistical vs practical importance.** With 569 samples, even small effects
  reach significance (symmetry error is an example). Effect sizes matter more
  than p-values here.
- **What the test shows.** Mann-Whitney compares the distributions of the two
  groups. It does not give a diagnostic rule: I did not train or validate a
  classifier, and nothing here supports clinical use.

## Conclusion

In this dataset, malignant breast mass samples have clearly larger and more
irregular nuclei than benign ones. Size and concavity features show the
strongest and most robust differences. The finding is an association from a
single, older, single-institution dataset and says nothing about causes or
clinical performance. The multiple-testing correction was not needed to separate
the strong signals from noise, but it keeps weak results like symmetry error
from being over-read.
