"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid contains a ring of atoms in its structure and retains characteristics of fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclic fatty acid, False otherwise
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

    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for hydrocarbon chains - adjusting the pattern to capture broader lengths
    # Altered to check for smaller sequences that could form longer chains through fragments
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No significant hydrocarbon chains found"

    return True, f"Contains {num_rings} ring structure(s) with fatty acid traits (carboxylic acid group and hydrocarbon chain)"