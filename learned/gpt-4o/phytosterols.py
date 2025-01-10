"""
Classifies: CHEBI:26125 phytosterols
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define a flexible sterol core pattern with stereochemistry consideration
    # A pattern for the steroid nucleus with possible stereochemistry variations
    sterol_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CC3C2CCC4(C3CCCC4)C')
    if not mol.HasSubstructMatch(sterol_pattern):
        return False, "No sterol core found"

    # Check for hydroxyl group typically attached to the sterol
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"

    # Count carbon atoms to validate side chain length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Phytosterols typically contain 27-30 carbons
    if c_count < 27 or c_count > 30:
        return False, f"Carbon count ({c_count}) not in typical phytosterol range"

    # Check for presence of side chain variations, could include double bonds
    side_chain_double_bond = Chem.MolFromSmarts('C=C')
    side_chain_matches = mol.GetSubstructMatches(side_chain_double_bond)
    
    if side_chain_matches:
        return True, "Contains sterol core with double bond(s) in side chain"

    return True, "Contains sterol core with typical side chains"

# Example usage:
# smiles = "[C@]1([C@@H](CCCC(C)C)C)([C@]2(CC[C@@]3([C@]4(CC[C@@H](C[C@]4(CC[C@]3([C@@]2(CC1)[H])[H])[H])O)C)[H])C)[H]"
# print(is_phytosterols(smiles))