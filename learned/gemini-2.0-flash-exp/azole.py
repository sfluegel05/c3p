"""
Classifies: CHEBI:68452 azole
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_azole(smiles: str):
    """
    Determines if a molecule is an azole based on its SMILES string.
    Azoles are monocyclic 5-membered heterocycles containing at least one nitrogen and
    possibly other heteroatoms (N, S, O).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an azole, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find 5-membered rings with at least 1 Nitrogen
    # This SMARTS pattern matches any 5-membered ring with at least one N and optionally other heteroatoms (O, S, N)
    azole_pattern = Chem.MolFromSmarts("[NX2,NX3;R5][*;R5]~[*;R5]~[*;R5]~[*;R5]") # 5-membered ring with at least one N
    
    ring_matches = mol.GetSubstructMatches(azole_pattern)

    if not ring_matches:
          return False, "No azole ring found"
    
    # Check the ring is monocyclic (no fused rings)
    for match in ring_matches:
      ring_atoms = [mol.GetAtomWithIdx(i) for i in match]
      if any(atom.IsInRing() for atom in ring_atoms): # Check if any atom in the ring is part of another ring.
        ring_systems = set()
        for atom in ring_atoms:
           for ring in atom.GetOwningRings():
             ring_systems.add(tuple(sorted(ring))) # add each owning ring
        if len(ring_systems) >1 : # if there is more than one ring in the system, it is not monocyclic.
          return False, "Found fused ring"


    return True, "Contains a monocyclic 5-membered heterocycle with at least one nitrogen."