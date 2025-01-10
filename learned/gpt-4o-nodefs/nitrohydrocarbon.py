"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon contains nitro groups (-NO2) attached to a hydrocarbon framework.

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

    # Define nitro group pattern
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)

    if len(nitro_matches) == 0:
        return False, "No nitro groups found"

    # Check carbon framework, ensuring majority are carbon atoms, common in hydrocarbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    total_atom_count = mol.GetNumAtoms()

    if c_count / total_atom_count < 0.4:
        return False, "Insufficient carbon atoms in framework"

    # Verify avoidance of extensive non-carbon ring/chain frameworks
    max_non_carbon = 0.2 * total_atom_count
    non_c_h_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in [1, 6, 7, 8])
    
    if non_c_h_count > max_non_carbon:
        return False, "Predominant non-carbon framework detected"

    return True, "Contains nitro groups attached to a hydrocarbon framework"