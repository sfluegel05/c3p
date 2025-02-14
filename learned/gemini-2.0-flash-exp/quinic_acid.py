"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid based on its SMILES string.
    A quinic acid is a cyclitol carboxylic acid, having a cyclohexane ring with multiple hydroxyl groups and a carboxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Core structure: cyclohexane with COOH
    # Use a SMARTS pattern that finds the core ring with a COOH group attached to a ring carbon.
    # The ring can have any combination of carbons, as long as there is one COOH attached and the ring is 6-membered
    quinic_acid_pattern = Chem.MolFromSmarts("[CX4]1([CX4][CX4][CX4][CX4][CX4]1)C(=O)[OX1H0]")
    if not mol.HasSubstructMatch(quinic_acid_pattern):
        return False, "Does not contain a cyclohexane ring with a carboxyl group directly attached"

    # Count the number of hydroxyl groups attached to the ring
    hydroxyl_count = 0
    ring_atoms = [atom for atom in mol.GetAtoms() if atom.IsInRing() and atom.GetAtomicNum() == 6]
    for ring_carbon in ring_atoms:
      for neighbor in ring_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
          hydroxyl_count += 1

    # Check for carboxylic acid group. Already verified in substructure match, but checking count just in case.
    carboxyl_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)[OX1H0]")))
    if carboxyl_count != 1:
      return False, "Must contain exactly one carboxylic acid group."
          
    if hydroxyl_count < 3:
      return False, f"Must contain at least three hydroxyl groups on the ring, found: {hydroxyl_count}"
    
    # Check that there is a 6-membered ring in molecule
    ring_info = mol.GetRingInfo()
    if not ring_info.IsAtomInRing(0):
         return False, "Molecule has not any ring."
    
    ring_sizes = ring_info.AtomRings()
    six_membered_ring = False
    for ring in ring_sizes:
      if len(ring) == 6:
        six_membered_ring = True
        break
    if not six_membered_ring:
        return False, "Molecule does not contain a 6-membered ring"

    return True, "Matches quinic acid pattern (cyclohexane ring with COOH and at least 3 hydroxyl groups)"