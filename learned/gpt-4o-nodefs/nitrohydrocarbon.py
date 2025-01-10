"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon contains nitro groups (-NO2) attached to a hydrocarbon framework.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for nitro groups [N+](=O)[O-]
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) == 0:
        return False, "No nitro groups found"

    # Check for predominant carbon and hydrogen framework
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)

    if c_count > 0 and h_count > 0:
        hydrocarbon_ratio = (c_count + h_count) / mol.GetNumAtoms()
        if hydrocarbon_ratio < 0.5:
            return False, "Insufficient hydrocarbon framework"
    else:
        return False, "No carbon or hydrogen atoms found"

    return True, "Contains nitro groups with a predominant hydrocarbon framework"