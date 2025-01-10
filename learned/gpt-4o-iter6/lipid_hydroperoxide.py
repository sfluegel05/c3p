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

    # Improved SMARTS for hydroperoxy groups, the [O][O][H] part may help
    hydroperoxy_pattern = Chem.MolFromSmarts("OO")
    if not mol.HasSubstructMatch(hydroperoxy_pattern):
        return False, "No hydroperoxy group found"

    # Detecting presence of long carbon chain regardless of exact carbon count
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long hydrocarbon chain detected"

    # Check for unsaturations which are typical in lipids
    unsaturation_pattern = Chem.MolFromSmarts("C=C")
    if mol.HasSubstructMatch(unsaturation_pattern):
        return True, "Contains hydroperoxy group(s) within unsaturated lipid-like chain"

    # Include ester groups detection as potential lipid characteristic
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(ester_pattern):
        return True, "Contains hydroperoxy group(s) within ester lipid-like structure"

    return False, "Missing typical structural features of a lipid despite presence of hydroperoxy group"