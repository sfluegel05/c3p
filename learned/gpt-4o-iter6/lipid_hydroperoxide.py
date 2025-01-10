"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy (-OOH) substituents.

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

    # Look for hydroperoxy groups (-OO[H])
    hydroperoxy_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2][OX2][OX2H]")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group found"

    # Assume lipid if it contains a long hydrocarbon chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Carbon chain too short for typical lipids: found {carbon_count} carbons"

    return True, "Contains hydroperoxy group(s) within a lipid structure"