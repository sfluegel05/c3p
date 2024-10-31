from rdkit import Chem
from rdkit.Chem import AllChem

def is_trihydroxyanthraquinone(smiles: str):
    """
    Determines if a molecule is a trihydroxyanthraquinone (anthraquinone with exactly 3 hydroxy groups)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a trihydroxyanthraquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check for anthraquinone core structure (two carbonyl groups connected by two benzene rings)
    anthraquinone_pattern = Chem.MolFromSmarts("O=C1c2ccccc2C(=O)c2ccccc12")
    if not mol.HasSubstructMatch(anthraquinone_pattern):
        return False, "No anthraquinone core structure found"
        
    # Count hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    if len(hydroxy_matches) != 3:
        return False, f"Found {len(hydroxy_matches)} hydroxy groups, need exactly 3"
        
    # Check if hydroxy groups are attached to aromatic carbons
    for hydroxy_match in hydroxy_matches:
        oxygen_atom = mol.GetAtomWithIdx(hydroxy_match[0])
        carbon_neighbors = [n for n in oxygen_atom.GetNeighbors() if n.GetSymbol() == 'C']
        if not carbon_neighbors or not carbon_neighbors[0].GetIsAromatic():
            return False, "One or more hydroxy groups not attached to aromatic ring"
            
    return True, "Valid trihydroxyanthraquinone structure found"
# Pr=1.0
# Recall=0.7857142857142857