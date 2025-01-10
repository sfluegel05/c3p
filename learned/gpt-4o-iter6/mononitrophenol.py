"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is defined as a nitrophenol where the phenol carries a
    single nitro substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        tuple: Reason for classification or potential compound details
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define phenol pattern
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1O")
    # Define nitro group pattern
    nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")
    
    # Check for phenol substructure
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol group found"
    
    # Check nitro groups in molecule
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) != 1:
        return False, f"Found {len(nitro_matches)} nitro groups, need exactly 1"
    
    # Ensure the nitro is part of the aromatic structure in phenol
    if not any(mol.GetAtomWithIdx(match[0]).GetIsAromatic() for match in nitro_matches):
        return False, "Nitro group is not attached to the aromatic ring of phenol"
    
    return True, "Molecule is a mononitrophenol with a single nitro group attached to the phenol ring"