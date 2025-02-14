"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is defined as a molecule with a cyclopenta[a]phenanthrene skeleton,
    partially or completely hydrogenated, with methyl groups at C-10 and C-13,
    and often an alkyl group at C-17.
    
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

    # Define the cyclopenta[a]phenanthrene core pattern using SMARTS (specific fused ring system)
    # The [C] in ring system are unspecified sp3 carbons
    # This pattern allows for some degree of unsaturation in the ring system, 
    # and is not very strict about the degree of saturation.
    steroid_core_smarts = "[C]1[C]2[C]3[C]([C]4[C]1[C]5[C]2[C]3[C]45)"
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)


    # Check if the molecule has the specific fused ring system
    if not mol.HasSubstructMatch(steroid_core_pattern):
       return False, "No cyclopenta[a]phenanthrene ring system found"

    # Find match of the core structure and use that to identify the relevant atoms.
    match = mol.GetSubstructMatch(steroid_core_pattern)
    atoms = [mol.GetAtomWithIdx(x) for x in match]

    # Now, let's define the carbon atoms where we expect methyl groups (C-10 and C-13, relative to the fused ring pattern)
    # Based on the SMARTS above, they are located in the following indices in the 'atoms' list. 
    # 0-1 is the start of 6 membered ring, 1-2 is the 2nd 6 membered ring.
    # 2-3 is the 3rd six membered ring, and then 3-4-5-0 is the 5 membered ring.
    # Indices 1 and 4 in this pattern would correspond to the carbons in which methyl groups are attached in steroids.
    # Note that this is based on the particular SMARTS we have defined.

    c10_idx = 1
    c13_idx = 4

    # Check for methyl groups at these positions
    methyl_count = 0
    for index, atom in enumerate(atoms):
        if index == c10_idx or index == c13_idx:
            for neighbor in atom.GetNeighbors():
              if neighbor.GetSymbol() == 'C' and len(neighbor.GetNeighbors()) == 1:
                  methyl_count += 1


    if methyl_count < 2:
       return False, f"Found only {methyl_count} methyl group(s) on positions 10 and 13, at least two are required."

    # Optional: Check for an alkyl group (at least one carbon) attached to C-17 (index 3 in the SMARTS match).
    # This is not strict requirement, but useful to add. We look for any chain of at least one carbon attached to index 3.
    alkyl_group_found = False
    c17_idx = 3
    for neighbor in atoms[c17_idx].GetNeighbors():
      if neighbor.GetSymbol() == 'C':
        alkyl_group_found = True
        break

    reason = "Contains the steroid tetracyclic ring system with methyl groups on C-10 and C-13."
    if alkyl_group_found:
        reason += " Also contains an alkyl group on C-17."

    return True, reason