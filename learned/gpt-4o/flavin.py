"""
Classifies: CHEBI:30527 flavin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_flavin(smiles: str):
    """
    Determines if a molecule is a flavin based on its SMILES string.
    A flavin is a derivative of dimethylisoalloxazine with a substituent on the 10th position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a flavin, False otherwise
        str: Reason for classification
    """
    # Convert SMILES string to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a looser SMARTS pattern for the dimethylisoalloxazine core
    core_pattern = Chem.MolFromSmarts("Cc1cc2nc3c([nH]c(=O)n(c3=O)c2cc1C)")
    
    # Check for core structure
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core flavin structure (dimethylisoalloxazine) missing"

    # Check for a substituent on the nitrogen atom (10th position equivalent)
    substituent_match = False
    core_matches = mol.GetSubstructMatches(core_pattern)
    
    for match in core_matches:
        substituted_nitrogen_idx = match[4]  # assuming reasonable index for substituent nitrogen
        # Check if any atom is connected to this nitrogen beyond core hydrogen
        atom = mol.GetAtomWithIdx(substituted_nitrogen_idx)
        for neighbor in atom.GetNeighbors():
            # If the neighbor is not hydrogen, it is a substituent
            if neighbor.GetAtomicNum() != 1:
                substituent_match = True
                break

    if not substituent_match:
        return False, "No substituent detected at the expected nitrogen position"
    
    return True, "Valid flavin structure detected"