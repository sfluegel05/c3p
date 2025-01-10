"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene is a C20 polyunsaturated fatty acid with four double bonds, including three conjugated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Leukotrienes typically have 20 carbon atoms
    if c_count != 20:
        return False, f"Incorrect number of carbon atoms: {c_count} (expected 20)"
    
    # Identify double bonds
    db_count = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    
    # Leukotrienes have four double bonds
    if db_count != 4:
        return False, f"Incorrect number of double bonds: {db_count} (expected 4)"
    
    # Check for conjugated double bonds (three should be conjugated)
    double_bond_pattern = Chem.MolFromSmarts('C=C')
    db_matches = mol.GetSubstructMatches(double_bond_pattern)
    
    # Verify at least three conjugated double bonds (simple connectivity check)
    conjugated_count = 0
    for i in range(len(db_matches) - 1):
        if db_matches[i][1] == db_matches[i + 1][0]:
            conjugated_count += 1
    
    if conjugated_count < 3:
        return False, "Less than 3 conjugated double bonds"

    return True, "Contains the characteristic C20 backbone with four double bonds, including three conjugated"