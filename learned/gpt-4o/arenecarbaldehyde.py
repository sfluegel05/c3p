"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is defined as an aldehyde in which the carbonyl group is attached to an aromatic moiety.

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
    
    # Check for presence of aromatic rings
    aromatic_atoms = [atom for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if not aromatic_atoms:
        return False, "No aromatic moiety found"
    
    # Look for aldehyde group pattern (-C(=O)H)
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not aldehyde_matches:
        return False, "No aldehyde group found"
    
    # Verify that the aldehyde is connected to an aromatic atom
    for match in aldehyde_matches:
        for idx in match:
            if mol.GetAtomWithIdx(idx).GetIsAromatic():
                return True, "Contains aldehyde group attached to an aromatic moiety"

    return False, "Aldehyde group not attached to aromatic moiety"