from rdkit import Chem
from rdkit.Chem import AllChem

def is_tetracyclic_antidepressant(smiles: str):
    """
    Determines if a molecule is a tetracyclic antidepressant based on structural features.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a tetracyclic antidepressant, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Get ring information
    rings = mol.GetRingInfo()
    ring_count = rings.NumRings()
    
    # Check for 4 rings
    if ring_count != 4:
        return False, f"Contains {ring_count} rings instead of required 4 rings"
        
    # Check for nitrogen atoms (typical in antidepressants)
    n_count = len([atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'N'])
    if n_count == 0:
        return False, "No nitrogen atoms found - required for antidepressant activity"
        
    # Get the rings as atom sets
    ring_sets = [set(ring) for ring in rings.AtomRings()]
    
    # Check for T-shaped arrangement by looking at ring connectivity
    # Two rings should share atoms with a central ring
    ring_connections = []
    for i in range(len(ring_sets)):
        for j in range(i+1, len(ring_sets)):
            if ring_sets[i].intersection(ring_sets[j]):
                ring_connections.append((i,j))
                
    # Count how many rings each ring is connected to
    connection_counts = {}
    for i,j in ring_connections:
        connection_counts[i] = connection_counts.get(i, 0) + 1
        connection_counts[j] = connection_counts.get(j, 0) + 1
        
    # For T-shape, one ring should connect to 3 others (central ring)
    if not any(count == 3 for count in connection_counts.values()):
        return False, "Rings not arranged in required T-shape configuration"
        
    # Check for aromatic rings (typical in tetracyclic antidepressants)
    aromatic_rings = 0
    for ring in ring_sets:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if all(atom.GetIsAromatic() for atom in atoms):
            aromatic_rings += 1
            
    if aromatic_rings < 2:
        return False, "Insufficient aromatic rings for tetracyclic antidepressant structure"

    # If all checks pass, classify as tetracyclic antidepressant
    return True, "Contains 4 rings in T-shape arrangement with nitrogen and aromatic systems"
# Pr=1.0
# Recall=0.6666666666666666