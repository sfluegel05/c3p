"""
Classifies: CHEBI:39362 mononitrophenol
"""
"""
Classifies: CHEBI:51088 mononitrophenol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol carrying a single nitro substituent at unspecified position.

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

    # Check for phenol group (aromatic ring with OH)
    phenol_pattern = Chem.MolFromSmarts("c1ccccc1[OH]")
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol group found"

    # Check for nitro groups ([N+]([O-])=O or N(=O)=O)
    nitro_pattern1 = Chem.MolFromSmarts("[N+](=[O-])=[O]")
    nitro_pattern2 = Chem.MolFromSmarts("[N](=O)=O")
    
    nitro_matches1 = len(mol.GetSubstructMatches(nitro_pattern1))
    nitro_matches2 = len(mol.GetSubstructMatches(nitro_pattern2))
    total_nitro = nitro_matches1 + nitro_matches2
    
    if total_nitro == 0:
        return False, "No nitro group found"
    elif total_nitro > 1:
        return False, f"Found {total_nitro} nitro groups, must have exactly one"

    # Check if the nitro group is attached to the phenol ring
    nitrophenol_pattern = Chem.MolFromSmarts("c1c([N+]([O-])=O)cccc1[OH]")
    nitrophenol_pattern2 = Chem.MolFromSmarts("c1c(N(=O)=O)cccc1[OH]")
    
    if not (mol.HasSubstructMatch(nitrophenol_pattern) or mol.HasSubstructMatch(nitrophenol_pattern2)):
        return False, "Nitro group not attached to phenol ring"

    return True, "Contains phenol ring with single nitro substituent"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:51088',
        'name': 'mononitrophenol',
        'definition': 'A nitrophenol that is phenol carrying a single nitro substituent at unspecified position.',
        'parents': ['CHEBI:33622']
    }
}