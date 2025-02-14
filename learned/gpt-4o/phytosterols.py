"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols found in plants, that have a similar core structure 
    to cholesterol, with variability in carbon side chains and the presence or 
    absence of double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a phytosterol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flexible sterol core pattern: a simple cyclopentanoperhydrophenanthrene nucleus
    sterol_core_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CC3C4CCC(C(C4)CC3)C2')
    if not mol.HasSubstructMatch(sterol_core_pattern):
        return False, "No sterol core found"
    
    # Check for hydroxyl group (optional for variants)
    hydroxyl_group = Chem.MolFromSmarts('[OX2H]')
    has_hydroxyl_group = mol.HasSubstructMatch(hydroxyl_group)
    
    # Check for variability in side chain - presence of a double bond
    side_chain_double_bond = Chem.MolFromSmarts('C=C')
    has_double_bond = mol.HasSubstructMatch(side_chain_double_bond)
    
    # Check overall carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (c_count >= 27 and c_count <= 35):
        return False, f"Carbon count ({c_count}) not in typical phytosterol range"
    
    # Justify the classification
    if has_hydroxyl_group and has_double_bond:
        return True, "Contains sterol core, hydroxyl group, and side chain variability with double bond"
    elif has_hydroxyl_group:
        return True, "Contains sterol core with hydroxyl group"
    elif has_double_bond:
        return True, "Contains sterol core with side chain variability (double bond)"
    
    return True, "Contains sterol core with typical side chain features"

# Example usage (for testing):
# smiles = "[C@]1([C@@H](CCCC(C)C)C)([C@]2(CC[C@@]3([C@]4(CC[C@@H](C[C@]4(CC[C@]3([C@@]2(CC1)[H])[H])[H])O)C)[H])C)[H]"
# print(is_phytosterols(smiles))