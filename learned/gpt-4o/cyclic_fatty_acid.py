"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

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
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxyl group (-COOH) or ester-like (-COO) patterns, allowing fatty acid derivatives
    carboxyl_or_ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H0,OX2H1]")
    if not mol.HasSubstructMatch(carboxyl_or_ester_pattern):
        return False, "No carboxyl or ester group found, not a recognized fatty acid derivative"

    # Calculate the number of carbon atoms to ensure a fatty acid characteristic
    c_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    # Ensure there's a reasonable length for what could be a fatty acid (usually > 6-8 carbons)
    if c_count < 8:
        return False, f"Too few carbon atoms ({c_count}), not a characteristic fatty acid or derivative"

    # Detect any cyclic structure
    if not mol.GetRingInfo().NumRings():
        return False, "No cyclic structure found"

    return True, "Contains cyclic structure with characteristics of a fatty acid or its derivative"