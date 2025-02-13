"""
Classifies: CHEBI:48030 tetrapeptide
"""
The previous program attempted to classify tetrapeptides by checking for the presence of 4 peptide bonds and 4 amino acid residues using SMARTS patterns. However, there are a few issues with this approach:

1. **False Positives**: The program only checks for the presence of 4 peptide bonds and 4 amino acid residues. This can lead to false positives for molecules that satisfy these conditions but are not actually tetrapeptides. For example, the program incorrectly classified several non-tetrapeptide molecules as tetrapeptides.

2. **Simplified Amino Acid Pattern**: The SMARTS pattern used to identify amino acid residues (`N[C@H](C)C(=O)`) is too simplistic and does not account for the diversity of amino acid side chains. This can lead to both false positives and false negatives.

3. **Molecular Weight Range**: The molecular weight range used (300-800 Da) may be too broad or narrow, leading to incorrect classifications.

4. **Elemental Composition**: While the elemental composition check is a good idea, it may not be sufficient on its own to distinguish tetrapeptides from other molecules.

To improve the program, we need to incorporate more specific structural features of tetrapeptides and potentially combine multiple criteria for a more robust classification.

Suggested improvements:

1. **Use a More Comprehensive Amino Acid Pattern**: Instead of a simplified pattern, use a more comprehensive SMARTS pattern or a set of patterns to accurately identify the diverse range of amino acid residues.

2. **Check for Specific Tetrapeptide Backbone**: Develop a SMARTS pattern or set of patterns to identify the specific backbone structure of a tetrapeptide, including the sequential arrangement of the four amino acid residues connected by peptide bonds.

3. **Incorporate Additional Structural Features**: Consider incorporating other structural features of tetrapeptides, such as the presence of specific functional groups (e.g., amino, carboxyl, amide), stereochemistry, or specific side chain patterns.

4. **Use a More Targeted Molecular Weight Range**: Refine the molecular weight range based on a larger dataset of known tetrapeptides to improve the accuracy of this criterion.

5. **Combine Multiple Criteria**: Instead of relying on a single criterion, combine multiple criteria (e.g., backbone pattern, amino acid residue patterns, molecular weight range, elemental composition) using logical operations to improve the overall classification accuracy.

6. **Consider Using Machine Learning Approaches**: For more complex classification tasks, consider using machine learning approaches that can learn the structural patterns of tetrapeptides from a large dataset of known examples.

By incorporating these improvements, the program should be able to more accurately classify tetrapeptides while minimizing false positives and false negatives.