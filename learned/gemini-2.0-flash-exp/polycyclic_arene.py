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

    # Get ring info
    ring_info = mol.GetRingInfo()
    
    # Get number of fused aromatic rings
    fused_rings_count = 0
    if ring_info.NumRings() > 1:
        for ring1_idx in range(ring_info.NumRings()):
           for ring2_idx in range(ring1_idx + 1, ring_info.NumRings()):
               ring1_atoms = set(ring_info.AtomRing(ring1_idx))
               ring2_atoms = set(ring_info.AtomRing(ring2_idx))
               if len(ring1_atoms.intersection(ring2_atoms)) > 0: # if rings share at least one atom
                   all_atoms_aromatic = True
                   for atom_idx in ring1_atoms.union(ring2_atoms):
                        atom = mol.GetAtomWithIdx(atom_idx)
                        if not atom.GetIsAromatic():
                            all_atoms_aromatic = False
                            break
                   if all_atoms_aromatic:
                       fused_rings_count+=1
        
    if fused_rings_count < 1:
        return False, "No fused aromatic rings detected"


    # Check for non-aromatic atoms other than hydrogens bonded to aromatic systems, and if those are too many or if they are not carbon
    non_aromatic_atoms = [atom for atom in mol.GetAtoms() if not atom.GetIsAromatic() and atom.GetAtomicNum() != 1]

    non_aromatic_atoms_bonded_to_aromatic = 0
    for atom in non_aromatic_atoms:
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIsAromatic():
                  non_aromatic_atoms_bonded_to_aromatic +=1
                  if atom.GetAtomicNum() != 6:
                      return False, f"Non-aromatic atom {atom.GetSymbol()} bonded to aromatic system."
                  break
    
    if len(non_aromatic_atoms) > 0 and non_aromatic_atoms_bonded_to_aromatic/len(non_aromatic_atoms) < 0.8:
         return False, f"Too many non-aromatic atoms, {len(non_aromatic_atoms)} found."
      
    
    # Check that at least 50% of the atoms are aromatic.
    total_atoms = mol.GetNumAtoms()
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())

    if total_atoms == 0:
        return False, "Molecule has no atoms"

    if aromatic_atoms/total_atoms < 0.5:
        return False, "Not enough aromatic atoms"

    return True, "Contains fused aromatic rings, classified as polycyclic arene"