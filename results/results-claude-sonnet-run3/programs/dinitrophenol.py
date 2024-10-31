from rdkit import Chem
from rdkit.Chem import AllChem

def is_dinitrophenol(smiles: str):
    """
    Determines if a molecule is a dinitrophenol (phenol with exactly 2 nitro groups).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a dinitrophenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phenol core
    phenol_pattern = Chem.MolFromSmarts('c1ccccc1O')
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "No phenol core found"
        
    # Count nitro groups
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    
    if len(nitro_matches) != 2:
        return False, f"Found {len(nitro_matches)} nitro groups, need exactly 2"
        
    # Get positions of nitro groups relative to OH
    oh_pattern = Chem.MolFromSmarts('cO')
    oh_match = mol.GetSubstructMatches(oh_pattern)[0]
    oh_carbon = oh_match[0]
    
    nitro_positions = []
    for match in nitro_matches:
        nitro_n = match[0]
        # Find carbon attached to nitro
        for bond in mol.GetAtomWithIdx(nitro_n).GetBonds():
            if bond.GetBeginAtomIdx() == nitro_n:
                other_atom = bond.GetEndAtomIdx()
            else:
                other_atom = bond.GetBeginAtomIdx()
            if mol.GetAtomWithIdx(other_atom).GetSymbol() == 'C':
                nitro_carbon = other_atom
                break
        
        # Get position relative to OH
        path = Chem.GetShortestPath(mol, oh_carbon, nitro_carbon)
        nitro_positions.append(len(path)-1)
        
    nitro_positions.sort()
    
    return True, f"Dinitrophenol with nitro groups at positions {nitro_positions}"
# Pr=1.0
# Recall=0.8888888888888888