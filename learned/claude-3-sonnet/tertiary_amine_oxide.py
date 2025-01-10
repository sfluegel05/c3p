"""
Classifies: CHEBI:134363 tertiary amine oxide
"""
"""
Classifies: tertiary amine N-oxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organic_group(mol, atom):
    """Helper function to check if an atom is part of an organic group"""
    if atom.GetAtomicNum() == 6:  # Carbon
        # Check if it's not bound to problematic atoms
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() not in [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]:  # H,C,N,O,F,P,S,Cl,Br,I
                return False
        return True
    return False

def is_tertiary_amine_oxide(smiles: str):
    """
    Determines if a molecule is a tertiary amine N-oxide based on its SMILES string.
    A tertiary amine N-oxide has a nitrogen atom bonded to three organic groups and 
    an oxygen atom with a formal charge separation (N+-O-).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary amine N-oxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for N+-O- pattern
    n_oxide_pattern = Chem.MolFromSmarts("[NX4+]-[O-]")
    if not mol.HasSubstructMatch(n_oxide_pattern):
        return False, "No N-oxide group found"
    
    # Get matches for N-oxide groups
    n_oxide_matches = mol.GetSubstructMatches(n_oxide_pattern)
    
    # Check each N-oxide nitrogen
    for match in n_oxide_matches:
        n_idx = match[0]  # Index of nitrogen atom
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check formal charge of N is +1
        if n_atom.GetFormalCharge() != 1:
            continue
            
        # Count organic neighbors
        organic_neighbors = 0
        non_oxide_neighbors = []
        
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetFormalCharge() == -1:
                continue  # Skip the oxide oxygen
            non_oxide_neighbors.append(neighbor)
            if is_organic_group(mol, neighbor):
                organic_neighbors += 1

        # Must have exactly 3 non-oxide neighbors
        if len(non_oxide_neighbors) != 3:
            continue
            
        # All three non-oxide neighbors must be organic groups
        if organic_neighbors != 3:
            continue
            
        # Check that nitrogen has no hydrogens
        if n_atom.GetTotalNumHs() != 0:
            continue
            
        # Additional structural checks
        
        # Exclude problematic patterns that shouldn't be classified as tertiary amine oxides
        problematic_patterns = [
            "[N+]1([O-])", # N-oxide in small ring
            "[N+]([O-])C=O", # N-oxide next to carbonyl
            "[N+]([O-])S", # N-oxide next to sulfur
            "[N+]([O-])P", # N-oxide next to phosphorus
        ]
        
        is_problematic = False
        for pattern in problematic_patterns:
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                is_problematic = True
                break
                
        if is_problematic:
            continue

        return True, "Found nitrogen with three organic groups and N-oxide"
    
    return False, "No nitrogen found with exactly three organic groups and N-oxide"