from rdkit import Chem
from rdkit.Chem import AllChem

def is_dimethylnaphthalene(smiles: str):
    """
    Determines if a molecule is a dimethylnaphthalene.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a dimethylnaphthalene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for naphthalene core
    # Should have exactly 10 carbons in fused aromatic rings
    rings = mol.GetRingInfo()
    
    # Find aromatic rings
    aromatic_rings = []
    for ring in rings.AtomRings():
        if len(ring) == 6:
            atoms = [mol.GetAtomWithIdx(i) for i in ring]
            if all(atom.GetIsAromatic() for atom in atoms):
                aromatic_rings.append(ring)
                
    if len(aromatic_rings) != 2:
        return False, "Not a naphthalene - must have exactly 2 aromatic rings"
        
    # Check if rings share atoms (are fused)
    ring1_atoms = set(aromatic_rings[0])
    ring2_atoms = set(aromatic_rings[1])
    if len(ring1_atoms.intersection(ring2_atoms)) != 2:
        return False, "Rings are not fused properly for naphthalene"
        
    # Get all ring atoms
    ring_atoms = ring1_atoms.union(ring2_atoms)
    
    # Check all ring atoms are carbon
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetSymbol() != 'C':
            return False, "Ring system contains non-carbon atoms"
            
    # Count methyl groups
    methyl_count = 0
    for atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in ring_atoms:
                # Check if it's a methyl group
                if (neighbor.GetSymbol() == 'C' and 
                    neighbor.GetDegree() == 1 and 
                    neighbor.GetTotalNumHs() == 3):
                    methyl_count += 1
                    
    if methyl_count != 2:
        return False, f"Found {methyl_count} methyl groups, must have exactly 2"
        
    positions = []
    for i, atom_idx in enumerate(sorted(ring_atoms)):
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in ring_atoms:
                if (neighbor.GetSymbol() == 'C' and 
                    neighbor.GetDegree() == 1 and 
                    neighbor.GetTotalNumHs() == 3):
                    positions.append(str(i+1))
    
    return True, f"{','.join(positions)}-dimethylnaphthalene"
# Pr=1.0
# Recall=1.0