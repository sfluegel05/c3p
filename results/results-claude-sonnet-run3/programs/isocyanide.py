from rdkit import Chem
from rdkit.Chem import AllChem

def is_isocyanide(smiles: str):
    """
    Determines if a molecule contains an isocyanide group (R-N≡C).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains isocyanide group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Look for [N+]#[C-] pattern
    isocyanide_pattern = Chem.MolFromSmarts('[N+]#[C-]')
    if mol.HasSubstructMatch(isocyanide_pattern):
        matches = mol.GetSubstructMatches(isocyanide_pattern)
        
        # Get atoms connected to nitrogen
        isocyanide_groups = []
        for match in matches:
            n_atom = mol.GetAtomWithIdx(match[0])
            neighbors = [x.GetSymbol() for x in n_atom.GetNeighbors() if x.GetIdx() != match[1]]
            
            if len(neighbors) == 1:  # N should have one neighbor besides C
                isocyanide_groups.append(neighbors[0])
        
        if isocyanide_groups:
            return True, f"Contains isocyanide group(s) with R = {', '.join(isocyanide_groups)}"
            
    return False, "No isocyanide group found"
# Pr=1.0
# Recall=1.0