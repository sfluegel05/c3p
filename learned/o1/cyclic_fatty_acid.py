"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
"""
Classifies: cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclic_fatty_acid(smiles: str):
    """
    Determines if a molecule is a cyclic fatty acid based on its SMILES string.
    A cyclic fatty acid is any fatty acid containing anywhere in its structure a ring of atoms.

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

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,C]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for ring structures
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No ring structures found"

    return True, f"Contains {num_rings} ring(s) and a carboxylic acid group"