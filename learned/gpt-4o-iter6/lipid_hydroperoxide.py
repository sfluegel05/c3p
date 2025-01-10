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
    
    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Modify to detect the presence of hydroperoxy group
    hydroperoxy_pattern = Chem.MolFromSmarts("[CX4,CX3,CX2]([OH1])OO")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group found"

    # Check for typical lipid-like structures: long carbon chains, ester groups, and unsaturations
    long_chain_pattern = Chem.MolFromSmarts("[C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0][C;R0]")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Look for ester functional groups typical in lipids
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    unsaturation_pattern = Chem.MolFromSmarts("C=C")
    
    if mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(unsaturation_pattern):
        return True, "Contains hydroperoxy group(s) within a lipid structure with unsaturation or ester group"

    return False, "Missing characteristics of a typical lipid structure despite presence of hydroperoxy group"