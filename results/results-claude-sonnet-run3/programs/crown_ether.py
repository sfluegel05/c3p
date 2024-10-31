from rdkit import Chem
from rdkit.Chem import AllChem

def is_crown_ether(smiles: str):
    """
    Determines if a molecule is a crown ether (cyclic compounds with repeating -O-CH2-CH2- units).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a crown ether, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule contains a ring
    rings = mol.GetRingInfo()
    if not rings.NumRings():
        return False, "No rings found - crown ethers must be cyclic"

    # Get the largest ring that contains oxygen
    oxygen_rings = []
    for ring in rings.AtomRings():
        oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetSymbol() == 'O')
        if oxygen_count > 0:
            oxygen_rings.append((ring, oxygen_count))
    
    if not oxygen_rings:
        return False, "No rings containing oxygen atoms found"
        
    largest_ring, num_oxygens = max(oxygen_rings, key=lambda x: len(x[0]))
    ring_atoms = set(largest_ring)

    # Check for proper -O-CH2-CH2- pattern
    for idx in largest_ring:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetSymbol() == 'O':
            # Check neighbors of oxygen
            neighbors = [n for n in atom.GetNeighbors() if n.GetIdx() in ring_atoms]
            if len(neighbors) != 2:
                return False, "Oxygen must have exactly two neighbors in ring"
            if not all(n.GetSymbol() == 'C' for n in neighbors):
                return False, "Oxygen must be connected to carbon atoms"
            
            # Check that neighboring carbons are sp3
            if any(n.GetIsAromatic() for n in neighbors):
                continue  # Skip this check for benzo-crown ethers
            if any(n.GetHybridization() != Chem.HybridizationType.SP3 for n in neighbors):
                return False, "Carbon atoms must be sp3 hybridized"

    # Count aromatic carbons for benzo-crown classification
    aromatic_carbons = sum(1 for idx in largest_ring 
                          if mol.GetAtomWithIdx(idx).GetIsAromatic())

    ring_size = len(largest_ring)
    
    # Validate minimum requirements
    if num_oxygens < 3:
        return False, "Crown ethers must have at least 3 oxygen atoms"
        
    if ring_size < 9:
        return False, "Ring size must be at least 9 for crown ethers"

    # Classify based on aromatic carbons
    if aromatic_carbons > 0:
        if aromatic_carbons == 6:
            return True, f"Benzo-{ring_size}-crown-{num_oxygens}"
        elif aromatic_carbons == 12:
            return True, f"Dibenzo-{ring_size}-crown-{num_oxygens}"
    
    # Regular crown ether
    return True, f"{ring_size}-crown-{num_oxygens}"
# Pr=1.0
# Recall=1.0