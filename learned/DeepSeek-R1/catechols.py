"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:33853 catechol
"""
from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol (contains an o-diphenol group) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains adjacent aromatic hydroxyl groups on a benzene ring, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings in the molecule
    rings = mol.GetRingInfo().AtomRings()
    
    for ring in rings:
        # Check if the ring is a benzene ring (6 atoms, all aromatic)
        if len(ring) != 6:
            continue
        
        is_aromatic = all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
        if not is_aromatic:
            continue
        
        # Convert ring to atom indices for easier adjacency checks
        ring_atoms = set(ring)
        
        # Check for adjacent hydroxyl groups
        for i, atom_idx in enumerate(ring):
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check current atom for hydroxyl (-OH)
            if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1:
                # Check next atom in the ring for adjacent hydroxyl
                next_idx = ring[(i+1) % 6]
                next_atom = mol.GetAtomWithIdx(next_idx)
                if next_atom.GetAtomicNum() == 8 and next_atom.GetTotalNumHs() >= 1:
                    return True, "Contains adjacent aromatic hydroxyl groups (o-diphenol) on benzene ring"
    
    # Check with SMARTS pattern as fallback for non-benzene aromatic systems
    catechol_pattern = Chem.MolFromSmarts("[OH]a-a[OH]")
    if mol.HasSubstructMatch(catechol_pattern):
        return True, "Contains adjacent aromatic hydroxyl groups (o-diphenol)"
    
    return False, "No o-diphenol component found"