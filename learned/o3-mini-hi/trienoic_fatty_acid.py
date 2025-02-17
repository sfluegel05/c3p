"""
Classifies: CHEBI:73155 trienoic fatty acid
"""
"""
Classifies: Trienoic fatty acid
Definition: Any polyunsaturated fatty acid that contains three double bonds.
"""

from rdkit import Chem

def is_trienoic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a trienoic fatty acid based on its SMILES string.
    A trienoic fatty acid is defined as a fatty acid (i.e. containing a carboxylic acid group)
    that contains exactly three carbon-carbon double bonds.
    
    Steps:
      1. Parse the input SMILES string.
      2. Check for the presence of a carboxylic acid group using a SMARTS pattern.
      3. Count the number of carbon-carbon double bonds (excluding the C=O from the carboxyl group).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a trienoic fatty acid, False otherwise.
        str: A reason describing the classification.
    """
    # Parse the SMILES string to create a molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the carboxylic acid group. We use the SMARTS "[CX3](=O)[OX2H1]",
    # which matches a carbonyl (C=O) bound to an -OH group.
    carboxylic_acid_smarts = "[CX3](=O)[OX2H1]"
    carboxylic_acid = Chem.MolFromSmarts(carboxylic_acid_smarts)
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "Missing carboxylic acid functional group (COOH)"
    
    # Count the number of double bonds between carbon atoms.
    # The double bonds in the carboxyl group (C=O) are not counted since one atom is oxygen.
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # Count the double bond only if both atoms are carbon.
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6:
                double_bond_count += 1

    # Verify that the molecule contains exactly three C=C double bonds.
    if double_bond_count != 3:
        return False, f"Found {double_bond_count} carbon-carbon double bonds; requires exactly 3 for trienoic fatty acid"

    return True, "Contains a carboxylic acid group and exactly 3 carbon-carbon double bonds (trienoic fatty acid phenotype)"