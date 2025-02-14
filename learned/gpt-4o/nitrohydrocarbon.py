"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more of the hydrogens has been replaced by nitro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is essentially a hydrocarbon
    # (only made up of C, H, and possibly attached nitro groups)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in {6, 1, 7, 8}:  # Carbon, Hydrogen, Nitrogen, Oxygen
            return False, f"Contains atom other than C, H, N, O: {atom.GetSymbol()}"

    # Check for nitro groups [N+](=O)[O-] attached to carbon atoms
    nitro_pattern = Chem.MolFromSmarts("[CX3][N+](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro group replacing hydrogen found"

    return True, "Contains hydrocarbon structure with one or more nitro groups replacing hydrogen"