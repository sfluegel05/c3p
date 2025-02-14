"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is defined as a molecule with a cyclopenta[a]phenanthrene skeleton,
    partially or completely hydrogenated, with methyl groups at C-10 and C-13.
    An alkyl group at C-17 is optional.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the cyclopenta[a]phenanthrene core pattern using SMARTS
    # The [R] in ring system are unspecified ring atoms
    steroid_core_smarts = "[R]1[R]2[R]3[R]([R]4[R]1[R]5[R]2[R]3[R]45)"
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)

    # Check if the molecule has the specific fused ring system
    if not mol.HasSubstructMatch(steroid_core_pattern):
       return False, "No cyclopenta[a]phenanthrene ring system found"

    # Check for at least 2 methyl groups attached to the ring system using SMARTS
    methyl_group_smarts = "[R][CH3]"
    methyl_group_pattern = Chem.MolFromSmarts(methyl_group_smarts)
    methyl_matches = mol.GetSubstructMatches(methyl_group_pattern)

    methyl_count = 0
    for match in methyl_matches:
      # Check that the ring atom that the methyl group is attached to
      #  is actually part of the steroid ring system.
      if mol.GetAtomWithIdx(match[0]).IsInRing():
          methyl_count += 1

    if methyl_count < 2:
      return False, f"Found only {methyl_count} methyl group(s) on the ring system, at least two are required."

    # Optional: Check for an alkyl group (at least one carbon) attached to the core.
    # Look for any chain of at least one carbon attached to any carbon of core
    alkyl_group_smarts = "[R](-[C])"
    alkyl_group_pattern = Chem.MolFromSmarts(alkyl_group_smarts)
    alkyl_matches = mol.GetSubstructMatches(alkyl_group_pattern)

    alkyl_group_found = False
    for match in alkyl_matches:
        if mol.GetAtomWithIdx(match[0]).IsInRing():
            alkyl_group_found = True
            break

    reason = "Contains the steroid tetracyclic ring system with at least two methyl groups on the core."
    if alkyl_group_found:
        reason += " Also contains an alkyl group on the core."

    return True, reason