"""
Classifies: CHEBI:48030 tetrapeptide
"""
Based on the previous code and the outcomes, it seems that the program is missing some key features to accurately classify tetrapeptides. Here are some potential issues and improvements:

1. **Molecular Weight Range**: The molecular weight range of 300-800 Da might be too broad or too narrow for tetrapeptides. Tetrapeptides typically have a molecular weight in the range of 400-800 Da, depending on the amino acid residues. Adjusting this range could improve accuracy.

2. **Backbone Pattern**: The SMARTS pattern used to identify the tetrapeptide backbone is too specific. It requires the amino acids to be connected in a specific order and does not account for variations in the backbone, such as cyclic peptides or branched peptides. A more flexible pattern or a combination of patterns might be necessary to capture all possible tetrapeptide structures.

3. **Amino Acid Residue Patterns**: The set of SMARTS patterns used to identify amino acid residues is comprehensive but might still miss some less common or modified residues. Additionally, the patterns do not account for protecting groups or other modifications that could be present in the SMILES strings.

4. **Stereochemistry**: The current program does not consider the stereochemistry of the amino acid residues, which is crucial for correctly identifying peptides. Tetrapeptides are typically composed of L-amino acids, so checking the stereochemistry could improve the classification accuracy.

5. **Handling Exceptional Cases**: Some of the false negatives in the outcomes might be due to exceptional cases or edge cases that the program does not handle correctly. For example, the program might struggle with peptides containing non-standard amino acids or unusual modifications.

Here are some potential improvements to address these issues:

- Refine the molecular weight range based on a more extensive analysis of tetrapeptide structures.
- Use a more flexible SMARTS pattern or a combination of patterns to identify the tetrapeptide backbone, accounting for variations such as cyclic peptides and branched peptides.
- Expand the set of SMARTS patterns for amino acid residues to include less common or modified residues, and consider handling protecting groups or other modifications.
- Incorporate stereochemistry checks to ensure that the identified amino acid residues have the correct stereochemistry (typically L-amino acids).
- Implement additional checks or exception handling for edge cases, such as peptides containing non-standard amino acids or unusual modifications.
- Consider using machine learning techniques or other more advanced methods to improve the classification accuracy, especially for complex or exceptional cases.

By addressing these issues and incorporating the suggested improvements, the program's ability to accurately classify tetrapeptides should be enhanced.