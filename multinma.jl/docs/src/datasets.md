# Bundled Datasets

multinma.jl ships with 9 datasets from published network meta-analyses.

## Listing and Loading

```julia
available_datasets()   # List all available dataset names
df = load_dataset("smoking")   # Load a specific dataset
```

## Dataset Descriptions

### `smoking` - Smoking Cessation (Hasselblad 1998)

24 studies, 4 treatments. Binary outcome (quit rates). A classic NMA dataset used in NICE TSD 4.

**Treatments:** No intervention, Self-help, Individual counselling, Group counselling

### `blocker` - Beta Blockers (Carlin 1992)

22 studies, 2 treatments. Binary outcome (mortality). A pairwise meta-analysis dataset.

**Treatments:** Control, Beta Blocker

### `thrombolytics` - Thrombolytic Treatments (Boland 2003)

50 studies, 9 treatments. Binary outcome (mortality after MI). A large NMA with complex network structure.

**Treatments:** SK, t-PA, Acc t-PA, SK + t-PA, r-PA, TNK, PTCA, UK, ASPAC

### `parkinsons` - Parkinson's Disease (TSD 2)

7 studies, 5 treatments. Continuous outcome. Available in both arm-based and contrast-based formats.

### `diabetes` - Diabetes (Elliott 2007)

22 studies, 6 treatments. Binary outcome. Demonstrates the complementary log-log (cloglog) link function.

**Treatments:** ACE Inhibitor, ARB, Beta Blocker, CCB, Diuretic, Placebo

### `statins` - Statin Cholesterol

19 studies. Binary outcome.

### `transfusion` - Blood Transfusion

6 studies. Binary outcome.

### `dietary_fat` - Dietary Fat Interventions

10 studies. Binary outcome.

### `atrial_fibrillation` - Atrial Fibrillation

26 studies. Mixed outcome types.

## Programmatic Datasets

Two datasets can also be created programmatically:

```julia
smoking = example_smoking()      # Smoking cessation (50 rows)
blocker = example_blocker()      # Beta blockers (10 rows, contrast-based)
```
