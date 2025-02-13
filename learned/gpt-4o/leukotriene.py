"""
Classifies: CHEBI:25029 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene is characterized as a C20 polyunsaturated fatty acid with specific conjugated double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check if the molecule has a minimum of 20 carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 20:
        return False, "Molecule does not have adequate icosanoid carbon chain length"
    
    # Check for the presence of four double bonds
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bonds) < 4:
        return False, f"Found {len(double_bonds)} double bonds, need at least 4"
    
    # Verify three consecutive conjugated double bonds (C=C-C=C-C=C)
    conjugated_bonds_pattern = Chem.MolFromSmarts("C=C-C=C-C=C")
    if not mol.HasSubstructMatch(conjugated_bonds_pattern):
        return False, "No conjugated double bonds pattern found"
    
    # Check for functional groups: hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_groups = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Assembly of the classification reason
    classification_reasons = [
        "Has adequate carbon chain length",
        "Has 4 double bonds, including conjugation pattern" if mol.HasSubstructMatch(conjugated_bonds_pattern) else "Double bond conjugation pattern missing",
        f"Contains {len(hydroxyl_groups)} hydroxyl groups" if hydroxyl_groups else "No hydroxyl groups detected"
    ]
    
    reason = "; ".join(filter(None, classification_reasons))
    
    # Returning the classification and reason
    return True, reason

# Testing example
# result, reason = is_leukotriene("CCCCCCCCC/C=C/C=C/C=C/C=C/CCCC(=O)O")
# print(result, reason)