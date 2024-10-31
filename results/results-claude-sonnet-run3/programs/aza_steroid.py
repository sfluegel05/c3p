from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetSymmSSSR

def is_aza_steroid(smiles: str):
    """
    Determines if a molecule is an aza-steroid (steroid with N replacing C in skeleton).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aza-steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic molecular properties
    num_atoms = mol.GetNumAtoms()
    if num_atoms < 20:  # Steroids typically have >20 atoms
        return False, "Too few atoms for steroid skeleton"
        
    # Check for presence of nitrogen
    nitrogen_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogen_atoms:
        return False, "No nitrogen atoms present"

    # Get all rings
    rings = list(GetSymmSSSR(mol))
    if len(rings) < 4:
        return False, "Does not contain minimum 4 rings required for steroid skeleton"

    # Find fused ring systems
    ring_systems = []
    processed = set()
    
    for i, ring1 in enumerate(rings):
        if i in processed:
            continue
            
        ring_system = set(ring1)
        processed.add(i)
        
        # Find all rings that share atoms with current ring system
        changed = True
        while changed:
            changed = False
            for j, ring2 in enumerate(rings):
                if j in processed:
                    continue
                if ring_system.intersection(ring2):
                    ring_system.update(ring2)
                    processed.add(j)
                    changed = True
                    
        ring_systems.append(ring_system)

    # Look for steroid ring system (largest fused ring system)
    if ring_systems:
        steroid_system = max(ring_systems, key=len)
        if len(steroid_system) < 17: # Typical steroid system has 17+ atoms
            return False, "Largest fused ring system too small for steroid"
    else:
        return False, "No fused ring systems found"

    # Check if any nitrogen is in the steroid ring system
    ring_nitrogens = [atom for atom in nitrogen_atoms if atom.GetIdx() in steroid_system]
    if ring_nitrogens:
        return True, f"Aza-steroid with {len(ring_nitrogens)} nitrogen(s) in steroid ring system"
    else:
        # Check if any nitrogen is substituted directly on the steroid ring system
        for atom_idx in steroid_system:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 7 and neighbor.GetIdx() not in steroid_system:
                    return True, "Aza-steroid with nitrogen substituent on steroid ring system"
                    
        return False, "No nitrogen atoms incorporated in steroid ring system"
# Pr=1.0
# Recall=0.8333333333333334