from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenanthridines(smiles: str):
    """
    Determines if a molecule is a phenanthridine (dibenzopyridine based on phenanthridine skeleton).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phenanthridine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Generate the aromatic ring information
    rings = mol.GetRingInfo()
    
    # Get all rings of size 6
    six_membered_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            ring_atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in ring_atoms):
                six_membered_rings.append(tuple(sorted(ring)))  # Use tuple instead of set
    
    if len(six_membered_rings) < 3:
        return False, "Does not contain required three fused aromatic rings"

    # Find nitrogen-containing rings (potential pyridine rings)
    n_containing_rings = []
    for ring in six_membered_rings:
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(atom.GetSymbol() == 'N' for atom in atoms):
            n_containing_rings.append(ring)

    if not n_containing_rings:
        return False, "No nitrogen-containing ring found"

    # Check for proper ring fusion pattern
    def are_rings_fused(ring1, ring2):
        return len(set(ring1).intersection(set(ring2))) >= 2

    def find_fused_systems(rings):
        systems = []
        used = set()
        
        for i, ring1 in enumerate(rings):
            if i in used:
                continue
                
            current_system = [i]
            changed = True
            
            while changed:
                changed = False
                for j, ring2 in enumerate(rings):
                    if j in current_system:
                        continue
                    if any(are_rings_fused(rings[k], ring2) for k in current_system):
                        current_system.append(j)
                        changed = True
            
            if len(current_system) >= 3:  # We need at least 3 fused rings
                systems.append([rings[i] for i in current_system])
            used.update(current_system)
            
        return systems

    fused_systems = find_fused_systems(six_membered_rings)
    
    if not fused_systems:
        return False, "No proper fusion pattern found"

    # Check if any fused system contains a nitrogen and proper connectivity
    for system in fused_systems:
        all_atoms = set()
        for ring in system:
            all_atoms.update(ring)
            
        n_atoms = [mol.GetAtomWithIdx(i) for i in all_atoms if mol.GetAtomWithIdx(i).GetSymbol() == 'N']
        
        if n_atoms and any(atom.GetIsAromatic() for atom in n_atoms):
            # Count carbons in the fused system
            c_atoms = len([mol.GetAtomWithIdx(i) for i in all_atoms if mol.GetAtomWithIdx(i).GetSymbol() == 'C'])
            
            if c_atoms >= 12:  # Phenanthridine core should have at least 12 carbons
                return True, "Contains phenanthridine core structure"

    return False, "Does not match phenanthridine structure requirements"
# Pr=1.0
# Recall=0.31746031746031744