"""
Classifies: CHEBI:37141 organobromine compound
"""
The previous program seems to have some issues that led to the low F1 score in classifying organobromine compounds:

1. **False Positives**: The program classifies some compounds as organobromine compounds even though they do not contain a carbon-bromine bond. This is likely because the program only checks for the presence of carbon and bromine atoms, and not specifically for a carbon-bromine bond. For example, compounds like "phloxine B" (SMILES: [Na+].[Na+].[O-]C(=O)c1c(Cl)c(Cl)c(Cl)c(Cl)c1-c1c2cc(Br)c([O-])c(Br)c2oc2c(Br)c(=O)c(Br)cc12) were wrongly classified as organobromine compounds, even though the bromine atoms are not directly bonded to carbon atoms.

2. **False Negatives**: The program fails to identify some compounds as organobromine compounds, even though they contain a carbon-bromine bond. This could be due to the way the program checks for carbon-bromine bonds. For example, compounds like "1-bromoindole" (SMILES: Brn1ccc2ccccc12) were missed, likely because the program checks for carbon-bromine bonds by iterating over the neighboring atoms of each bromine atom, and in this case, the bromine atom is directly bonded to a nitrogen atom, not a carbon atom.

To improve the program, we can consider the following modifications:

1. **Use SMARTS patterns**: Instead of iterating over atoms and checking their neighbors, we can use SMARTS patterns to match specific substructures in the molecule. SMARTS patterns can be used to directly identify carbon-bromine bonds.

2. **Handle aromaticity**: Some compounds, like "1-bromoindole," have aromatic rings where the carbon-bromine bond is not explicitly represented in the SMILES string. We need to handle such cases by considering aromaticity and bond types.

3. **Consider bond orders**: The current program does not consider bond orders, which could lead to false positives or false negatives in certain cases. For example, a compound with a triple bond between carbon and bromine should not be considered an organobromine compound.

4. **Utilize RDKit's built-in functionality**: RDKit provides various functions and descriptors that can be used to identify specific substructures and bond types, which could simplify the implementation and improve accuracy.

By addressing these issues, we can enhance the program's performance in classifying organobromine compounds accurately.