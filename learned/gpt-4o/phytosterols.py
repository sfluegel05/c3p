"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol with variations mainly 
    in carbon side chains and sometimes presence of a double bond.

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

    # Sterol core: cyclopentanoperhydrophenanthrene nucleus with possible stereochemistry
    sterol_pattern = Chem.MolFromSmarts('[C@@H]1CC[C@H]2[C@H](C1)[C@@H]3CC[C@@]4(C)[C@@H](CC[C@]4(C)[C@@]3(C)[C@@H]2C)O')
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol core found"

    # Check for presence of hydroxyl group, characteristic in sterols
    hydroxyl_group = Chem.MolFromSmarts('[CX4][OX2H]')
    if not mol.HasSubstructMatch(hydroxyl_group):
        return False, "No hydroxyl group found"

    # Check for side chain variability
    side_chain_double_bond = Chem.MolFromSmarts('CCC=C')
    if mol.HasSubstructMatch(side_chain_double_bond):
        return True, "Contains sterol core with double bond(s) in side chain"

    # Check general carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 27 or c_count > 35:
        return False, f"Carbon count ({c_count}) not in typical phytosterol range"

    return True, "Contains sterol core with typical side chain features"

# Example usage (for testing):
# smiles = "[C@]1([C@@H](CCCC(C)C)C)([C@]2(CC[C@@]3([C@]4(CC[C@@H](C[C@]4(CC[C@]3([C@@]2(CC1)[H])[H])[H])O)C)[H])C)[H]"
# print(is_phytosterols(smiles))