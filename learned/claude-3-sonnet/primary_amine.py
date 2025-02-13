"""
Classifies: CHEBI:32877 primary amine
"""
Analyzing the previous code and the outcomes, it appears that the program correctly identifies most primary amines, but also wrongly classifies some molecules as primary amines, leading to false positives. Additionally, it misses some primary amines, resulting in false negatives.

The main issue with the previous code is that it relies solely on the presence of a nitrogen atom with one substituent (excluding hydrogen) to classify a molecule as a primary amine. This approach is too simplistic and fails to account for other structural features and bonding patterns that may or may not be consistent with a primary amine.

To improve the program, we need to incorporate additional checks and criteria to better distinguish primary amines from other nitrogen-containing compounds. Here are some suggestions:

1. **Check for aromatic nitrogen atoms**: Primary amines typically have an aliphatic (non-aromatic) nitrogen atom. We should exclude molecules with aromatic nitrogen atoms, as these are more likely to be aromatic amines or heterocycles rather than primary amines.

2. **Check for multiple nitrogen atoms**: Primary amines should have only one nitrogen atom. Molecules with multiple nitrogen atoms are likely to be other types of compounds, such as amines with other substituents or heterocycles.

3. **Check for double or triple bonds to nitrogen**: Primary amines have a single bond between the nitrogen atom and its substituent. Molecules with double or triple bonds to the nitrogen atom should be excluded, as these are more likely to be imines, nitriles, or other nitrogen-containing functional groups.

4. **Check for other functional groups**: Primary amines should not have certain functional groups, such as carbonyl groups (e.g., amides, esters, or carboxylic acids), which would indicate a different type of compound.

5. **Check for specific structural patterns**: Primary amines often have a characteristic structural pattern, such as a alkyl chain or a benzyl group attached to the nitrogen atom. We could use substructure matching to identify these patterns and increase the confidence in the classification.

6. **Consider molecular weight and atom counts**: Primary amines typically have a relatively low molecular weight and a specific range of carbon, nitrogen, and hydrogen atom counts. We could set thresholds or ranges for these properties to further refine the classification.

By incorporating these additional checks and criteria, we can significantly improve the accuracy of the primary amine classification program and reduce the number of false positives and false negatives.