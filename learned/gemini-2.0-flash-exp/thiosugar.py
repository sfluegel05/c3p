"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar is a carbohydrate derivative in which one or more of the oxygens
    or hydroxy groups of the parent carbohydrate is replaced by sulfur or -SR,
    where R can be hydrogen or any group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of sulfur atoms
    sulfur_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16]
    if not sulfur_atoms:
      return False, "No sulfur atoms found"
    
    # Check for a carbohydrate-like core (ring with multiple hydroxyls or oxygens)
    # Relaxed SMARTS pattern to detect 5 or 6-membered rings with carbons and oxygens
    sugar_ring_pattern = Chem.MolFromSmarts("[C,O;R5,R6]1[C,O;R5,R6][C,O;R5,R6][C,O;R5,R6][C,O;R5,R6]1")
    sugar_rings_matches = mol.GetSubstructMatches(sugar_ring_pattern)
    if not sugar_rings_matches:
        return False, "No carbohydrate-like ring structure detected"
    
    # Check for at least 3 oxygen atoms in the ring or attached to it
    oxygens_in_or_attached_to_ring = 0
    for ring_match in sugar_rings_matches:
      for atom_idx in ring_match:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 8:
          oxygens_in_or_attached_to_ring += 1
        for neighbor in atom.GetNeighbors():
           if neighbor.GetAtomicNum() == 8:
               oxygens_in_or_attached_to_ring +=1
    
    if oxygens_in_or_attached_to_ring < 3:
        return False, "Too few oxygens in or attached to ring"

    # Check if sulfur is connected to a sugar ring, or is part of SR group.
    sulfur_connected_to_ring = False
    for s_atom in sulfur_atoms:
        for neighbor in s_atom.GetNeighbors():
            for ring_match in sugar_rings_matches:
               if neighbor.GetIdx() in ring_match:
                   sulfur_connected_to_ring = True
                   break
            if sulfur_connected_to_ring:
              break
        if sulfur_connected_to_ring:
            break
    if not sulfur_connected_to_ring:
      return False, "Sulfur not directly attached to a sugar ring."

    #Check for multiple sugar rings with sulfur connections
    if len(sugar_rings_matches) > 1 and sulfur_connected_to_ring:
         return True, "Multiple carbohydrate rings with sulfur present"

    return True, "Sulfur is attached to a carbohydrate-like ring."