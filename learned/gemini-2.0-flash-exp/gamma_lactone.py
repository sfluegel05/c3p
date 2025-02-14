"""
Classifies: CHEBI:37581 gamma-lactone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is a lactone with a five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the gamma-lactone ring: C(=O)-O-C-C
    gamma_lactone_pattern = Chem.MolFromSmarts('[CX3](=[OX1])-[OX2]-[CX4]-[CX4]-')

    # Check if the molecule has the gamma-lactone pattern
    matches = mol.GetSubstructMatches(gamma_lactone_pattern)
    if not matches:
          return False, "No gamma-lactone ring pattern found."
    
    for match in matches:
      # Check the size of the ring
      oxygen_atom_index = match[1]
      carbon1_atom_index = match[2]
      carbon2_atom_index = match[3]
      oxygen_atom = mol.GetAtomWithIdx(oxygen_atom_index)

      # Get the list of neighbor atoms
      neighbors = [atom.GetIdx() for atom in oxygen_atom.GetNeighbors()]
      if carbon1_atom_index not in neighbors:
          return False, "Incorrect gamma-lactone ring pattern: connectivity issue."
      if carbon2_atom_index not in neighbors:
          return False, "Incorrect gamma-lactone ring pattern: connectivity issue."

      # Check if the neighbours of the carbonyl carbon also form a 5-membered ring with the oxygen atoms.
      carbonyl_carbon_index = match[0]
      carbonyl_carbon_atom = mol.GetAtomWithIdx(carbonyl_carbon_index)
      carbonyl_carbon_neighbours = [atom.GetIdx() for atom in carbonyl_carbon_atom.GetNeighbors()]
      if oxygen_atom_index not in carbonyl_carbon_neighbours:
         return False, "Incorrect gamma-lactone ring pattern: connectivity issue."
      
      num_carbons_in_ring = 0
      for idx in match:
          if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6:
              num_carbons_in_ring += 1
      if num_carbons_in_ring != 4:
          return False, "Incorrect ring size. Gamma-lactones must have four carbon atoms in the ring"
         
    return True, "Molecule contains a gamma-lactone ring."