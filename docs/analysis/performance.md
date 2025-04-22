# Performance

## Macro stats

![img](comparison/model_comparison_macro.png)

## Micro stats

![img](comparison/model_comparison_micro.png)

## Rank distrubution across models

For each model, we ranked each class by f1, and then plotted the
distrubution of ranks across models. This gives a sense of which
classes are consistently learnable, and classes for which models
varied in how well they learned them

![img](comparison/rank_distribution_top_50.png)


## Generalization from training to test data

![img](ensemble/scatterplot-train_f1-f1.png)

## Similarity of models

To explore the relative contribution of each approach, we examined the
correlation between f1 scores on each class for each pair of methods
(supplementary methods). In general, there was high correlation
between different C3POs, but less correlation between C3POs and
Chebifier. This indicates that deep learning and program learning
approaches are likely complementary:

![img](comparison/scatter_matrix.png)
