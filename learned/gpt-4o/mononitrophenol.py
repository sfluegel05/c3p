"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol with a single nitro group attached directly on the same benzene ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for phenol with a single nitro group
    # Match phenol group (c-O) with a single nitro group attached anywhere
    mononitrophenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")

    # Find substructure matches
    phenol_matches = mol.GetSubstructMatches(mononitrophenol_pattern)
    
    # Find all nitro groups and ensure there is only one
    nitro_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'N' and mol.HasSubstructMatch(nitro_pattern, atom.GetIdx()))

    # Check conditions
    if len(phenol_matches) > 0 and nitro_count == 1:
        return True, "Contains phenol group with a single nitro group on the same benzene ring"
    if len(phenol_matches) <= 0:
        return False, "No phenol group found"
    if nitro_count != 1:
        return False, f"Invalid number of nitro groups: found {nitro_count}"

    return False, "Does not meet mononitrophenol criteria"