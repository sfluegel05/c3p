"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
from rdkit import Chem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon in which one or more hydrogens have been replaced by nitro groups.
    
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

    # Define nitro group pattern -[N+](=O)[O-]
    nitro_pattern = Chem.MolFromSmarts("[NX3](=O)[O-]")
    if not mol.HasSubstructMatch(nitro_pattern):
        return False, "No nitro group found"

    # Check carbon backbone by verifying that the primary atoms are Carbon and possibly Hydrogen
    non_hydrocarbon_atoms = [atom.GetAtomicNum() for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6, 7, 8)]
    if non_hydrocarbon_atoms:
        return False, "Contains non-hydrocarbon atoms"

    return True, "Molecule is a nitrohydrocarbon with one or more nitro groups attached to a carbon backbone"