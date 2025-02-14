"""
Classifies: CHEBI:33848 polycyclic arene
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene based on its SMILES string.
    A polycyclic arene consists of fused aromatic rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polycyclic arene, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of aromatic rings
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic rings found"
    
    # Check if there are multiple aromatic rings
    aromatic_count = len(mol.GetSubstructMatches(aromatic_ring_pattern))
    if aromatic_count < 2:
        return False, "Less than 2 aromatic rings found, not a polycyclic arene"

    # Identify all rings
    ring_atoms = []
    for atom in mol.GetAtoms():
        if atom.IsInRing() and atom.GetIsAromatic():
            ring_atoms.append(atom.GetIdx())
    
    if not ring_atoms:
      return False, "No aromatic rings detected"

    # Group atoms into rings using the ring info.
    rings = []
    
    for bond in mol.GetBonds():
        if bond.IsInRing():
          ring_indices = set()
          for atom in [bond.GetBeginAtom().GetIdx(), bond.GetEndAtom().GetIdx()]:
              if atom not in ring_atoms:
                  continue
              
              found = False
              for ring_idx in range(len(rings)):
                  if atom in rings[ring_idx]:
                    found = True
                    break
              if not found:
                ring_indices.add(atom)
                
          if len(ring_indices) > 0:
              
            new_ring = list(ring_indices)
            new_ring_set = set(new_ring)
            
            
            
            # check for fusion: if at least one atom is already in other ring, fuse the rings
            merged = False
            for ring_idx in range(len(rings)):
                if new_ring_set.intersection(set(rings[ring_idx])):
                    rings[ring_idx].extend(new_ring)
                    rings[ring_idx] = list(set(rings[ring_idx])) # Remove duplicates.
                    merged=True
                    break
            if not merged:
                rings.append(new_ring) # if no existing ring is found, create a new ring.
                
    
    # Count fused aromatic rings
    fused_rings_count = 0
    if len(rings) > 1:
        for ring1_idx in range(len(rings)):
           for ring2_idx in range(ring1_idx + 1, len(rings)):
               ring1_atoms = set(rings[ring1_idx])
               ring2_atoms = set(rings[ring2_idx])
               if len(ring1_atoms.intersection(ring2_atoms)) > 0: # if rings share at least one atom
                  fused_rings_count += 1
        
    if fused_rings_count < 1:
      return False, "No fused aromatic rings detected"

    # Check for non-aromatic atoms (excluding hydrogens)
    non_aromatic_atoms = [atom for atom in mol.GetAtoms() if not atom.GetIsAromatic() and atom.GetAtomicNum() != 1]
    
    if len(non_aromatic_atoms) > mol.GetNumAtoms() * 0.3:
         return False, f"Too many non-aromatic atoms, {len(non_aromatic_atoms)} found."

    # Check that at least 50% of the atoms are aromatic.
    total_atoms = mol.GetNumAtoms()
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())

    if total_atoms == 0:
        return False, "Molecule has no atoms"

    if aromatic_atoms/total_atoms < 0.5:
        return False, "Not enough aromatic atoms"


    return True, "Contains fused aromatic rings, classified as polycyclic arene"