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

    # Look for carboxylic acid group or similar structures
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    ester_pattern = Chem.MolFromSmarts("C(=O)O")  # to include ester-like patterns
    if not mol.HasSubstructMatch(carboxylic_acid_pattern) and not mol.HasSubstructMatch(ester_pattern):
        return False, "No carboxylic acid or ester group found"

    # Check for hydrocarbon chains, allowing for flexibility in length and interruptions
    hydrocarbon_chain_pattern = Chem.MolFromSmarts("[CH2]~[CH,CH2]*~[CH2]")
    if not mol.HasSubstructMatch(hydrocarbon_chain_pattern):
        return False, "No significant hydrocarbon chains found"

    # Consider inclusion of exocyclic functionalities or specific structural features if appropriate
    # For example, epoxy, hydroxyl, or double bonds adjacent to rings

    return True, f"Contains {num_rings} ring structure(s) with fatty acid traits (carboxylic acid or ester group and hydrocarbon chain)"