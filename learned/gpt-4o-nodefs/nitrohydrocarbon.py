"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon contains nitro groups (-NO2) attached to a primarily hydrocarbon framework.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the nitro group pattern
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)

    if len(nitro_matches) == 0:
        return False, "No nitro groups found"

    # Count carbon and other critical atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    total_atom_count = mol.GetNumAtoms()

    # Allow isotope variations but keep focus on regular atoms
    isotope_count = sum(1 for atom in mol.GetAtoms() if atom.GetIsotope() != 0)

    # Calculate adjusted values removing isotopes
    effective_carbon_ratio = (c_count / (total_atom_count - isotope_count))

    if effective_carbon_ratio < 0.45:
        return False, "Insufficient carbon atoms in framework"

    # Ensure no excessive functional groups besides nitro (Focus mainly on hydrocarbon character)
    non_hydrocarbon_counts = total_atom_count - c_count - n_count - o_count - isotope_count
    if non_hydrocarbon_counts > 2:
        return False, "Contains heteroatoms or functional groups not typical in a hydrocarbon framework"

    return True, "Contains nitro groups attached to a primarily hydrocarbon framework"