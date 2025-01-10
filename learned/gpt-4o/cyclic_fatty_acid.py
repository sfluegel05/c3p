"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid features at least one ring structure and retains characteristics of fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a cyclic fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify ring structures
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No ring structure found"

    # Look for carboxylic acid or ester group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    contains_acid_or_ester = mol.HasSubstructMatch(carboxylic_acid_pattern) or mol.HasSubstructMatch(ester_pattern)
    if not contains_acid_or_ester:
        return False, "No carboxylic acid or ester group found"

    # Check for hydrocarbon chains with more flexibility
    hydrocarbon_pattern1 = Chem.MolFromSmarts("[CH2]~[CH,CH2]*~[CH2]")  # Linear or interrupted chains
    hydrocarbon_pattern2 = Chem.MolFromSmarts("C~[CH2,CH,OH]")  # Account for common functional groups adjacent
    has_hydrocarbon_chain = mol.HasSubstructMatch(hydrocarbon_pattern1) or mol.HasSubstructMatch(hydrocarbon_pattern2)
    if not has_hydrocarbon_chain:
        return False, "No significant hydrocarbon chains found"

    return True, f"Contains {num_rings} ring structure(s) with fatty acid traits (carboxylic acid or ester group and hydrocarbon chain)"