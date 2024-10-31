from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiocyanates(smiles: str):
    """
    Determines if a molecule is a thiocyanate (RSC#N).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a thiocyanate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Look for SC#N substructure
    thiocyanate_pattern = Chem.MolFromSmarts('[S;X2][C;X2]#[N;X1]')
    matches = mol.GetSubstructMatches(thiocyanate_pattern)
    
    if not matches:
        return False, "No thiocyanate (SC#N) group found"
        
    # Check that sulfur is connected to carbon (R group)
    for match in matches:
        sulfur_idx = match[0]
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)
        
        # Get neighbors of sulfur excluding the cyano carbon
        neighbors = [n for n in sulfur_atom.GetNeighbors() if n.GetIdx() != match[1]]
        
        if len(neighbors) != 1:
            return False, "Sulfur must have exactly one R group attached"
            
        r_group = neighbors[0]
        if r_group.GetSymbol() != 'C':
            return False, "R group must be carbon-based"
            
    num_matches = len(matches)
    if num_matches == 1:
        return True, "Contains one thiocyanate (RSC#N) group"
    else:
        return True, f"Contains {num_matches} thiocyanate (RSC#N) groups"
# Pr=1.0
# Recall=1.0