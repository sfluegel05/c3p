"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid should contain features with a cyclic component and fatty acid-like moieties,
    which traditionally include long carbon chains and terminal carboxylic acid groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure there are ring structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 1:
        return False, "No ring structures found"

    # Look for carboxylic acid groups
    carboxylic_acid_group = Chem.MolFromSmarts("C(=O)O")
    has_carboxylic_acid = mol.HasSubstructMatch(carboxylic_acid_group)

    # Detect common cyclic structures such as furans and epoxides
    fragment_patterns = [
        Chem.MolFromSmarts("c1ccoc1"),  # Furan
        Chem.MolFromSmarts("[O]1[C][C][O]1"),  # Epoxide or similar cyclic ethers
    ]

    has_alternative_cycle = any(mol.HasSubstructMatch(pattern) for pattern in fragment_patterns)
    
    if has_carboxylic_acid and (ring_info.NumRings() > 1 or has_alternative_cycle):
        return True, "Contains both cyclic and fatty acid features, classifying as a cyclic fatty acid"
    
    return False, "Does not fully fit the cyclic fatty acid criteria based on initial checks"

# Example usage:
# result, reason = is_cyclic_fatty_acid("OC(=O)CCCC[C@H]1CCC=C1")
# print(result, reason)