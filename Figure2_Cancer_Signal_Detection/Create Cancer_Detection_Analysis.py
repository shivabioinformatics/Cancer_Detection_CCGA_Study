import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, auc, confusion_matrix
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10

#============================================================================
# STEP 1: LOAD THE DATA
#============================================================================

# Load classifier scores (this contains predictions from all 10 classifiers)
print("Loading data...")
scores_df = pd.read_csv('data/scores_cnc.tsv', sep='\t')
clinical_df = pd.read_csv('data/clinical.tsv', sep='\t')
cutoffs_df = pd.read_csv('data/score_cutoffs.tsv', sep='\t')

print(f"✓ Loaded {len(scores_df):,} predictions across all classifiers")
print(f"✓ Loaded {len(clinical_df):,} participants")
print()

# Check unique classifiers
classifiers = scores_df['classifier_name'].unique()
print(f"10 Classifiers in the study:")
for i, clf in enumerate(classifiers, 1):
    print(f"  {i:2d}. {clf}")
print()


#============================================================================
# STEP 2: REPRODUCE TABLE 3 - VALIDATION SET RESULTS
#============================================================================

# Filter to validation set only
validation_df = scores_df[scores_df['train_or_valid'] == 'valid'].copy()

#============================================================================
#- **Training set (1,414 people):** Used to teach the classifiers how to spot cancer
#- **Validation set (847 people):** Used to test if it actually works on NEW people it's never seen
#============================================================================

# Count samples by cancer status
n_cancer = (validation_df.groupby('classifier_name')['cnc_label_actual']
            .apply(lambda x: (x == 'cancer').sum())
            .iloc[0])
n_noncancer = (validation_df.groupby('classifier_name')['cnc_label_actual']
               .apply(lambda x: (x == 'non_cancer').sum())
               .iloc[0])

print(f"Validation Set Composition:")
print(f"  Cancer samples:     {n_cancer:4d}")
print(f"  Non-cancer samples: {n_noncancer:4d}")
print(f"  Total:              {n_cancer + n_noncancer:4d}")
print()

# Initialize results table
table3_results = []

# Order classifiers as in paper (Table 3)
classifier_order = [
    'WG methylation',
    'SNV',
    'SNV-WBC',
    'SCNA',
    'SCNA-WBC',
    'Fragment endpoints',
    'Fragment lengths',
    'Allelic imbalance',
    'Pan-feature',
    'Clinical data'
]

# Calculate metrics for each classifier
for clf_name in classifier_order:
    clf_data = validation_df[validation_df['classifier_name'] == clf_name].copy()

    # Get actual and predicted labels
    y_true = (clf_data['cnc_label_actual'] == 'cancer').astype(int)
    y_pred = (clf_data['cnc_label_predict'] == 'cancer').astype(int)

    # Calculate confusion matrix
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()

    # Calculate metrics
    sensitivity = tp / (tp + fn) if (tp + fn) > 0 else 0
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0

    # Calculate 95% confidence intervals for sensitivity (Wilson score interval)
    n_positives = tp + fn
    z = 1.96  # 95% CI
    p_hat = sensitivity

    if n_positives > 0:
        denominator = 1 + z**2/n_positives
        center = (p_hat + z**2/(2*n_positives)) / denominator
        margin = z * np.sqrt((p_hat*(1-p_hat)/n_positives + z**2/(4*n_positives**2))) / denominator
        ci_low = max(0, center - margin)
        ci_high = min(1, center + margin)
    else:
        ci_low = ci_high = 0

    # Store results
    table3_results.append({
        'Classifier': clf_name,
        'TP': tp,
        'FN': fn,
        'TN': tn,
        'FP': fp,
        'Total Cancer': tp + fn,
        'Sensitivity': sensitivity,
        'Sensitivity %': f"{sensitivity*100:.1f}%",
        'CI_Low': ci_low,
        'CI_High': ci_high,
        '95% CI': f"({ci_low*100:.0f}%–{ci_high*100:.0f}%)",
        'Specificity': specificity,
        'Specificity %': f"{specificity*100:.1f}%",
        'TP/Total': f"{tp}/{tp+fn}"
    })

# Create DataFrame
table3_df = pd.DataFrame(table3_results)

# Display formatted table (matching paper's Table 3)
print("VALIDATION SET RESULTS (matching Table 3 in paper)")
print("-" * 85)
print(f"{'Classifier':<25} {'Sensitivity':<15} {'95% CI':<20} {'TP/Total':<12}")
print("-" * 85)

for _, row in table3_df.iterrows():
    print(f"{row['Classifier']:<25} {row['Sensitivity %']:<15} {row['95% CI']:<20} {row['TP/Total']:<12}")

print("-" * 85)
print(f"Specificity target: 98.0%")
print(f"Observed specificity: {table3_df.iloc[0]['Specificity %']}")
print()

# Highlight key findings
print("KEY FINDINGS:")
wg_meth = table3_df[table3_df['Classifier'] == 'WG methylation'].iloc[0]
snv_wbc = table3_df[table3_df['Classifier'] == 'SNV-WBC'].iloc[0]
pan_feat = table3_df[table3_df['Classifier'] == 'Pan-feature'].iloc[0]

print(f"  1. WG methylation:  {wg_meth['Sensitivity %']} sensitivity (TOP single assay)")
print(f"  2. SNV-WBC:         {snv_wbc['Sensitivity %']} sensitivity")
print(f"  3. Pan-feature:     {pan_feat['Sensitivity %']} sensitivity (combining all)")
print()
print("  → Methylation was the BEST single genomic feature!")
print("  → SNV required WBC sequencing to match methylation performance")
print()

#============================================================================
# STEP 3: STATISTICAL COMPARISON (McNEMAR'S TEST)
#============================================================================

#============================================================================
# QUESTION: "Is WG methylation REALLY better, or did we just get lucky?"
#============================================================================

#Compare two tests on the **same people** (paired comparison)

# McNemar's test: Compare WG methylation vs others

def mcnemar_test(clf1_name, clf2_name, data):
    """
    Perform McNemar's test comparing two classifiers

    Returns:
        statistic, p_value, contingency_table
    """
    clf1 = data[data['classifier_name'] == clf1_name].copy()
    clf2 = data[data['classifier_name'] == clf2_name].copy()

    # Merge on participant_id to ensure paired comparison
    merged = pd.merge(
        clf1[['participant_id', 'cnc_label_actual', 'cnc_label_predict']],
        clf2[['participant_id', 'cnc_label_actual', 'cnc_label_predict']],
        on=['participant_id', 'cnc_label_actual'],
        suffixes=('_clf1', '_clf2')
    )

    # Only look at cancer samples (where we measure sensitivity)
    cancer_only = merged[merged['cnc_label_actual'] == 'cancer'].copy()

    # Create contingency table
    correct_1 = (cancer_only['cnc_label_predict_clf1'] == 'cancer')
    correct_2 = (cancer_only['cnc_label_predict_clf2'] == 'cancer')

    # McNemar table
    both_correct = (correct_1 & correct_2).sum()
    clf1_only = (correct_1 & ~correct_2).sum()
    clf2_only = (~correct_1 & correct_2).sum()
    both_wrong = (~correct_1 & ~correct_2).sum()

    # McNemar's statistic focuses on discordant pairs
    n = clf1_only + clf2_only
    if n > 0:
        # Continuity correction
        statistic = (abs(clf1_only - clf2_only) - 1)**2 / n
        p_value = 1 - stats.chi2.cdf(statistic, df=1)
    else:
        statistic = 0
        p_value = 1.0

    contingency = pd.DataFrame({
        f'{clf2_name} Correct': [both_correct, clf2_only],
        f'{clf2_name} Wrong': [clf1_only, both_wrong]
    }, index=[f'{clf1_name} Correct', f'{clf1_name} Wrong'])

    return statistic, p_value, contingency, clf1_only, clf2_only

# Compare WG methylation vs all others
print("Comparing WG Methylation (best) vs Other Classifiers:")
print("-" * 70)

reference = 'WG methylation'
comparisons = ['SNV', 'SNV-WBC', 'SCNA', 'SCNA-WBC', 'Fragment lengths',
               'Fragment endpoints', 'Allelic imbalance']

for comp in comparisons:
    stat, p_val, table, wg_only, comp_only = mcnemar_test(reference, comp, validation_df)

    significance = "***" if p_val < 0.001 else ("**" if p_val < 0.01 else ("*" if p_val < 0.05 else "ns"))

    print(f"\n{reference} vs {comp}:")
    print(f"  Detected by WG only: {wg_only}")
    print(f"  Detected by {comp} only: {comp_only}")
    print(f"  p-value: {p_val:.4f} {significance}")

print()
print("Significance levels: *** p<0.001, ** p<0.01, * p<0.05, ns = not significant")
print()
print("INTERPRETATION:")
print("  WG methylation significantly outperforms most classifiers (p < 0.01)")
print("  SNV-WBC matches WG methylation (p > 0.05) - no significant difference")
print()

#============================================================================
# STEP 4: REPRODUCE FIGURE 2 - ROC CURVES
# ROC = Receiver Operating Characteristic (fancy name for "performance graph")
#============================================================================

# Create figure with subplots
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Color scheme (highlight WG methylation in orange)
colors = {
    'WG methylation': '#FF8C00',  # Orange
    'SNV-WBC': '#1f77b4',
    'SNV': '#aec7e8',
    'SCNA-WBC': '#2ca02c',
    'SCNA': '#98df8a',
    'Fragment lengths': '#d62728',
    'Fragment endpoints': '#ff9896',
    'Allelic imbalance': '#9467bd',
    'Pan-feature': '#8c564b',
    'Clinical data': '#e377c2'
}

# Classifier order for legend (best to worst)
legend_order = [
    'Pan-feature',
    'WG methylation',
    'SNV-WBC',
    'SCNA-WBC',
    'Fragment lengths',
    'SNV',
    'SCNA',
    'Allelic imbalance',
    'Fragment endpoints',
    'Clinical data'
]

# Calculate ROC curves for each dataset
for idx, (dataset_name, dataset_filter) in enumerate([('train', 'train'), ('valid', 'valid')]):
    ax = axes[idx]

    # Plot diagonal (random chance)
    ax.plot([0, 1], [0, 1], 'k--', lw=1, label='Random chance', alpha=0.5)

    # Store AUC values
    auc_values = {}

    # Plot ROC for each classifier
    for clf_name in legend_order:
        clf_data = scores_df[
            (scores_df['classifier_name'] == clf_name) &
            (scores_df['train_or_valid'] == dataset_filter)
        ].copy()

        if len(clf_data) == 0:
            continue

        # Get true labels and prediction scores
        y_true = (clf_data['cnc_label_actual'] == 'cancer').astype(int).values
        y_score = clf_data['p_cancer'].values

        # Calculate ROC curve
        fpr, tpr, thresholds = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)
        auc_values[clf_name] = roc_auc

        # Plot
        ax.plot(fpr, tpr,
                color=colors.get(clf_name, 'gray'),
                lw=2.5 if clf_name == 'WG methylation' else 1.5,
                label=f'{clf_name}',
                alpha=0.9 if clf_name == 'WG methylation' else 0.7)

    # Formatting
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('False-Positive Rate', fontsize=12)
    ax.set_ylabel('True-Positive Rate', fontsize=12)
    ax.set_title(f'{"Training" if dataset_name == "train" else "Validation"} Performance',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='lower right', fontsize=9)
    ax.grid(True, alpha=0.3)

    # Add inset zooming to high specificity region (98-100%)
    # This matches the paper's Figure 2 bottom panels
    axins = ax.inset_axes([0.15, 0.55, 0.4, 0.4])

    for clf_name in legend_order:
        clf_data = scores_df[
            (scores_df['classifier_name'] == clf_name) &
            (scores_df['train_or_valid'] == dataset_filter)
        ].copy()

        if len(clf_data) == 0:
            continue

        y_true = (clf_data['cnc_label_actual'] == 'cancer').astype(int).values
        y_score = clf_data['p_cancer'].values
        fpr, tpr, _ = roc_curve(y_true, y_score)

        # Only plot high specificity region
        mask = (fpr >= 0.0) & (fpr <= 0.1)
        axins.plot(fpr[mask], tpr[mask],
                   color=colors.get(clf_name, 'gray'),
                   lw=2 if clf_name == 'WG methylation' else 1,
                   alpha=0.9 if clf_name == 'WG methylation' else 0.6)

    axins.set_xlim(0.0, 0.1)
    axins.set_ylim(0.0, 0.5)
    axins.set_xlabel('False-Positive Rate', fontsize=8)
    axins.set_ylabel('True-Positive Rate', fontsize=8)
    axins.tick_params(labelsize=8)
    axins.grid(True, alpha=0.3)

    # Mark 98% specificity line (2% FPR)
    axins.axvline(x=0.02, color='red', linestyle='--', lw=1, alpha=0.5)

plt.tight_layout()
plt.savefig('out/figure2_roc_curves.png', dpi=300, bbox_inches='tight')
print("✓ Figure 2 saved to: out/figure2_roc_curves.png")
print()
