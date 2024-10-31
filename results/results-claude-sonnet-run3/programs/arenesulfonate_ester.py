from rdkit import Chem
from rdkit.Chem import AllChem

def is_arenesulfonate_ester(smiles: str):
    """
    Determines if a molecule is an arenesulfonate ester.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an arenesulfonate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Look for S(=O)(=O)O pattern
    sulfonyl_pattern = Chem.MolFromSmarts('S(=O)(=O)O')
    if not mol.HasSubstructMatch(sulfonyl_pattern):
        return False, "No sulfonate group found"

    # Get matches for sulfonate group
    matches = mol.GetSubstructMatches(sulfonyl_pattern)
    
    for match in matches:
        sulfur_idx = match[0]
        oxygen_idx = match[3]
        
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        
        # Check if sulfur is connected to an aromatic ring
        aromatic_ring = False
        for neighbor in sulfur_atom.GetNeighbors():
            if neighbor.GetIsAromatic():
                aromatic_ring = True
                break
                
        if not aromatic_ring:
            continue
            
        # Check if oxygen is connected to carbon (ester)
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetSymbol() == 'C' and neighbor.GetIdx() != sulfur_idx:
                # Found arenesulfonate ester pattern
                return True, "Contains arenesulfonate ester group (Ar-SO2-O-R)"
                
    return False, "No arenesulfonate ester group found"
# Pr=1.0
# Recall=0.8461538461538461