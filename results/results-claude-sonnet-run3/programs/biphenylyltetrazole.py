from rdkit import Chem
from rdkit.Chem import AllChem

def is_biphenylyltetrazole(smiles: str):
    """
    Determines if a molecule contains a biphenyl ring system substituted by a tetrazole ring.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a biphenylyltetrazole, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # SMARTS patterns
    biphenyl_pattern = Chem.MolFromSmarts('c1ccccc1-c1ccccc1') # Biphenyl core
    tetrazole_pattern = Chem.MolFromSmarts('c1nnn[nH]1')  # Tetrazole ring
    
    # Check for biphenyl
    if not mol.HasSubstructMatch(biphenyl_pattern):
        return False, "No biphenyl substructure found"
        
    # Check for tetrazole
    if not mol.HasSubstructMatch(tetrazole_pattern):
        return False, "No tetrazole substructure found"
        
    # Get matches
    biphenyl_matches = mol.GetSubstructMatches(biphenyl_pattern)
    tetrazole_matches = mol.GetSubstructMatches(tetrazole_pattern)
    
    # Check if tetrazole is connected to biphenyl
    biphenyl_atoms = set()
    for match in biphenyl_matches:
        biphenyl_atoms.update(match)
        
    for tetrazole_match in tetrazole_matches:
        # Get neighbors of tetrazole atoms
        tetrazole_neighbors = set()
        for atom_idx in tetrazole_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                tetrazole_neighbors.add(neighbor.GetIdx())
                
        # Check if any tetrazole neighbor is part of biphenyl
        if tetrazole_neighbors & biphenyl_atoms:
            return True, "Contains biphenyl substituted with tetrazole"
            
    return False, "Tetrazole not connected to biphenyl system"
# Pr=1.0
# Recall=1.0