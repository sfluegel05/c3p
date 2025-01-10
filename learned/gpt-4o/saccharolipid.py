"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is characterized by the presence of both a lipid and a carbohydrate moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carbohydrate moiety: Presence of sugar-like structure (e.g., cyclic ethers)
    carbohydrate_pattern = Chem.MolFromSmarts('[C;R](O)[C;R](O)')  # Simple representation for cyclic sugars
    if not mol.HasSubstructMatch(carbohydrate_pattern):
        return False, "No carbohydrate moiety found"

    # Look for lipid moiety: Long carbon chains typically ending in a carboxyl group
    lipid_pattern = Chem.MolFromSmarts('C(=O)[O][C@H](C)CCCCCCCC')  # Simple representation for lipids
    if not mol.HasSubstructMatch(lipid_pattern):
        return False, "No lipid moiety found"
    
    # Confirm linkage between the carbohydrate and lipid regions
    linkage_pattern = Chem.MolFromSmarts('[C@H](O)[C;R](O)C(=O)[O][C@H](C)CCCCCCCC')  # Example linkage pattern
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No linkage between carbohydrate and lipid found"

    return True, "Contains linked carbohydrate and lipid moieties, consistent with a saccharolipid"

# Example usage:
# result, reason = is_saccharolipid("input_smiles_here")
# print(result, reason)