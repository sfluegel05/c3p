"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenolic compound with an aromatic ring carrying a single nitro group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a mononitrophenol, False otherwise
        tuple: Reason for classification or error
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # SMARTS pattern to detect a phenol ring (aromatic ring with OH group)
    phenol_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")
    
    # SMARTS pattern for a nitro group attached to an aromatic carbon
    nitro_pattern = Chem.MolFromSmarts("c[N+](=O)[O-]")
    
    # Check if the molecule contains a phenol pattern
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenolic structure found"
    
    # Find all nitro groups
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if len(nitro_matches) == 0:
        return False, "No nitro group found"
    
    # Verify only one nitro group and directly linked to a phenol ring
    for match in nitro_matches:
        connected_aromatic = any([atom.GetSymbol() == 'O' and atom.IsInRing() for atom in mol.GetBondBetweenAtoms(match[0], match[1]).GetAtoms()])
        if connected_aromatic:
            return True, "Valid mononitrophenol with appropriately attached nitro group"
    
    return False, "The molecule has improper nitro linkage or multiple nitro groups"