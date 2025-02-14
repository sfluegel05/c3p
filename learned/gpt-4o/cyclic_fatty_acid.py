"""
Classifies: CHEBI:59238 cyclic fatty acid
"""
from rdkit import Chem

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

    # Check for carboxyl group (-COOH or COO)
    carboxyl_pattern = Chem.MolFromSmarts("[C](=O)[O;H1,-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found, not a fatty acid"

    # Calculate the number of carbons
    c_count = sum(atom.GetAtomicNum() == 6 for atom in mol.GetAtoms())
    if c_count < 8:
        return False, f"Too few carbon atoms ({c_count}), not a fatty acid"

    # Detect any cyclic structure
    if not mol.GetRingInfo().NumRings():
        return False, "No cyclic structure found"

    # Check for density of non-fatty-acid related heteroatoms
    possible_impurities = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 8]) / (c_count + 1)
    if possible_impurities > 0.2:  # arbitrary threshold for significant 'non-fatty' impurities
        return False, "Excessive non-fatty-acid related heteroatoms"

    return True, "Contains cyclic structure with characteristic fatty acid components"