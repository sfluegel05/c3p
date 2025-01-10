"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is identified by a benzene ring with one hydroxyl group and one nitro group.

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

    # Look for benzene ring with a hydroxyl group
    phenol_pattern = Chem.MolFromSmarts("Oc1ccccc1")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol group found"

    # Look for a nitro group (in RDKit, often represented as -[#7](=[O])-[O])
    nitro_group_pattern = Chem.MolFromSmarts("[$([NX3](=O)=O)]")
    nitro_matches = mol.GetSubstructMatches(nitro_group_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, need exactly 1"

    return True, "Contains a phenol group with one nitro group attached"