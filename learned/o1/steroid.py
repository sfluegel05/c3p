"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: Steroids
"""

from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is defined as a compound based on the cyclopenta[a]phenanthrene skeleton,
    which consists of three fused six-membered rings and one five-membered ring fused together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a steroid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid core SMARTS pattern
    steroid_core_smarts = """
    [#6]1([#6][#6][#6]2[#6][#6][#6]3[#6][#6][#6]([#6]3[#6][#6][#6]4[#6][#6][#6]1[#6][#6][#6]24)[#6][#6][#6]5[#6][#6][#6]([#6][#6][#6]5)[#6][#6]2)
    """
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Failed to define steroid core SMARTS pattern"

    # Check for the steroid core in the molecule
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Steroid core skeleton not found in molecule"

    # Optionally, check for methyl groups attached to ring carbons
    methyl_count = 0
    for atom in mol.GetAtoms():
        if atom.IsInRing() and atom.GetAtomicNum() == 6:
            # Iterate over neighbors to find methyl groups
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1:
                    # Neighbor is a methyl group (CH3)
                    methyl_count += 1

    if methyl_count < 2:
        # Steroids usually have methyl groups at C-10 and C-13
        return True, f"Steroid core found, but only {methyl_count} methyl groups attached to ring carbons"

    return True, "Steroid core skeleton found with methyl groups attached to ring carbons"