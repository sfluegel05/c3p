"""
Classifies: CHEBI:86315 methyl sulfide
"""
from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    A methyl sulfide is an aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the methyl sulfide SMARTS pattern
    methyl_sulfide_pattern = Chem.MolFromSmarts("[SX2]-[CX4H3]")

    # Check for substructure match of a methyl sulfide
    if not mol.HasSubstructMatch(methyl_sulfide_pattern):
        return False, "Molecule does not contain a methyl sulfide group"

    # Iterate through all the sulfur atoms, and check that each one matches the required criteria
    for sulfur in mol.GetSubstructMatches(Chem.MolFromSmarts("[S]")):
      sulfur_atom = mol.GetAtomWithIdx(sulfur[0])
      if sulfur_atom.GetExplicitValence() != 2:
          continue
      
      # Get neighbors of the sulfur atom
      neighbors = [atom for atom in sulfur_atom.GetNeighbors()]
      if len(neighbors) != 2:
            continue
      
      carbon_atoms = [atom for atom in neighbors if atom.GetAtomicNum() == 6]
      if len(carbon_atoms) != 2:
          continue

      # Check that at least one is a methyl
      methyl_match = False
      for carbon in carbon_atoms:
        if carbon.GetExplicitValence() == 4 and carbon.GetTotalNumHs() == 3:
              methyl_match = True
              break
      if not methyl_match:
        continue
      
      #Check that the carbon atoms are aliphatic
      aliphatic_match = True
      for carbon in carbon_atoms:
          if carbon.GetIsAromatic():
              aliphatic_match = False
              break
      if aliphatic_match:
          return True, "Molecule is an aliphatic methyl sulfide"

    return False, "Molecule does not contain an aliphatic methyl sulfide"