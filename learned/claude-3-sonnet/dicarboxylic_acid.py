"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: CHEBI:35652 dicarboxylic acid
Any carboxylic acid containing two carboxy groups.
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid contains at least two carboxyl groups (-C(=O)O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxyl group pattern (-C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H0]")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)

    # Check if there are at least 2 carboxyl groups
    if len(carboxyl_matches) < 2:
        return False, f"Found {len(carboxyl_matches)} carboxyl groups, need at least 2"

    # Additional checks for common dicarboxylic acid substructures
    # Check for cyclic or aromatic systems with carboxyl groups
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        ring_carboxyl_count = sum(sum(mol.GetAtomWithIdx(idx).IsInRingSize(6) for idx in match) for match in carboxyl_matches)
        if ring_carboxyl_count >= 2:
            return True, "Contains at least two carboxyl groups in a cyclic or aromatic system"

    # Check for specific functional group patterns (e.g., alpha-keto acids, amino acids)
    # ...

    return True, "Contains at least two carboxyl groups (-C(=O)O)"