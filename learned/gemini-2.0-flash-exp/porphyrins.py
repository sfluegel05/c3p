"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: CHEBI:26389 porphyrin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is defined by a macrocycle consisting of four pyrrole rings linked by methine bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for pyrrole rings and methine bridges (more general this time)
    pyrrole_pattern = Chem.MolFromSmarts("[n;H0]1cc[c]c1") # pyrrole
    methine_pattern = Chem.MolFromSmarts("[CX3]") # methine

    # Find pyrrole rings
    pyrrole_matches = mol.GetSubstructMatches(pyrrole_pattern)
    if len(pyrrole_matches) != 4:
        return False, f"Found {len(pyrrole_matches)} pyrrole rings, require 4."
    
    # Find methine groups
    methine_matches = mol.GetSubstructMatches(methine_pattern)
    if len(methine_matches) < 4:
        return False, f"Found {len(methine_matches)} methine groups, require at least 4."

    # Extract atom indices from the matches
    pyrrole_atoms = [match[0] for match in pyrrole_matches] # take only the first atom of each pyrrole match
    methine_atoms = [match[0] for match in methine_matches] # take only the first atom of each methine match


    # Check connectivity: each pyrrole must be connected to two methines.
    # For this, I need to check if any pyrrole atom is connected to the found methine atoms
    # In a porphyrin core, each pyrrole is connected by a methine bridge, forming a macrocycle.
    # First, find all the bonds that connect a pyrrole and a methine
    
    connecting_bonds = 0 # count how many connecting bonds are present
    for p_atom in pyrrole_atoms:
      for neighbor in mol.GetAtomNeighbors(p_atom):
          if neighbor.GetIdx() in methine_atoms:
            connecting_bonds +=1

    if connecting_bonds < 8: # need 2 bonds per pyrrole, 4 pyrroles, then 4*2=8
        return False, "Pyrrole rings and methine bridges are not correctly connected"
      
    return True, "Contains four pyrrole rings linked by methine bridges."