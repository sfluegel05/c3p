"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde using its SMILES string.
    An arenecarbaldehyde is an aldehyde where the carbonyl group is attached to an aromatic moiety.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Detect aromaticity in the molecule
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    
    # Check for aldehyde group pattern (-C(=O)H)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    # Check if carbon of aldehyde is directly attached to an aromatic atom
    for match in aldehyde_matches:
        aldehyde_carbon_idx = match[0]
        aldehyde_carbon_atom = mol.GetAtomWithIdx(aldehyde_carbon_idx)
        
        for neighbor in aldehyde_carbon_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                aromatic_system = False
                for aromatic_neighbor in neighbor.GetNeighbors():
                    if aromatic_neighbor.GetIsAromatic():
                        aromatic_system=True
                if not aromatic_system:
                    continue
                return True, "Contains aldehyde group directly attached to an aromatic moiety"
    
    return False, "Aldehyde group not attached to aromatic moiety"