"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is a lipid molecule with one or more hydroperoxy (-OOH) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroperoxy group pattern (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("[OX2H][OX1]")  # Improved SMARTS to match -OOH correctly
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group (-OOH) found"

    # Ensure hydroperoxy group is connected correctly to a carbon chain (generic lipid structure)
    for match in mol.GetSubstructMatches(hydroperoxy_pattern):
        atom_index = match[1]  # Oxygen in -OOH
        neighbors = mol.GetAtomWithIdx(atom_index).GetNeighbors()
        # Confirm one of the neighbors is carbon, and part of a long chain
        if any(neighbor.GetAtomicNum() == 6 for neighbor in neighbors):
            return True, "Contains hydroperoxy groups attached to a lipid backbone"

    return False, "Hydroperoxy group not attached correctly to carbon backbone"