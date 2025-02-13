"""
Classifies: CHEBI:25676 oligopeptide
"""
The previous program attempts to classify molecules as oligopeptides based on several criteria, such as the number of amino acid residues, the presence of peptide bonds and amino acid sidechains, and the molecular weight range. However, there are several issues with the current implementation that could lead to the observed low performance.

1. **Amino acid residue count**: The SMARTS pattern used to count amino acid residues ([N;X3;H2,H1&!$(N(-C=O)-O-C=O)][C;X4][C;X3](=O)[N;X3]) may not cover all possible variations of amino acid residues. For example, it might miss residues with modified or unusual side chains, or residues with unusual bond orders or charges.

2. **Peptide bond and sidechain detection**: The SMARTS patterns used to detect peptide bonds and amino acid sidechains are very general and may match unintended substructures. More specific patterns might be needed to accurately identify these features.

3. **Molecular weight range**: The molecular weight range used (300-2000 Da) may be too restrictive or not accurate for all oligopeptides. Oligopeptides can have a wider range of molecular weights depending on their sequence and modifications.

4. **Lack of context**: The current program does not consider the context of the peptide within the larger molecule. For example, it might incorrectly classify larger peptides or proteins as oligopeptides if they contain a short peptide fragment.

To improve the classification performance, the following steps could be taken:

1. **Improve amino acid residue detection**: Use a more comprehensive set of SMARTS patterns or a machine learning model trained on known amino acid residues to more accurately identify and count them.

2. **Refine peptide bond and sidechain detection**: Use more specific SMARTS patterns or machine learning models to detect peptide bonds and amino acid sidechains, taking into account different modifications and unusual structures.

3. **Adjust molecular weight range**: Analyze a larger dataset of known oligopeptides to determine a more accurate molecular weight range, or consider using machine learning to learn the appropriate range from the data.

4. **Consider context**: Incorporate logic to identify the peptide fragment within the larger molecule and consider its context. This could involve detecting the presence of other functional groups or substructures that are incompatible with oligopeptides.

5. **Use machine learning**: Consider training a machine learning model on a large dataset of known oligopeptides and non-oligopeptides, using a combination of molecular descriptors and structural features as input. This could potentially learn more complex patterns and relationships than a rule-based approach.

6. **Combine with other methods**: Integrate the current approach with other methods, such as similarity searches or substructure matching against known oligopeptides, to improve the overall classification accuracy.

It's worth noting that the task of accurately classifying oligopeptides based solely on their SMILES strings is challenging, and a combination of different approaches might be necessary to achieve high performance.