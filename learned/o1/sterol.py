"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: CHEBI:15889 sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3β-hydroxy steroid whose skeleton is closely related to cholestan-3-ol,
    allowing for additional carbon atoms in the side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid backbone SMARTS pattern (cyclopentanoperhydrophenanthrene ring system)
    steroid_pattern = Chem.MolFromSmarts('C1CC2CCC1CCC3CCC4(C)CCCC4C3C2')  # Simplified steroid core

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found (cyclopentanoperhydrophenanthrene)"

    # Find the matches for the steroid backbone
    steroid_matches = mol.GetSubstructMatches(steroid_pattern)
    if not steroid_matches:
        return False, "No steroid backbone matches found"

    # Check for 3β-hydroxyl group
    # In steroids, the 3-position is on ring A
    hydroxyl_pattern = Chem.MolFromSmarts('[C;R1]([O;H1])[C;R1]')  # Carbon with hydroxyl attached in a ring

    # Search for hydroxyl group attached to ring carbon in the steroid backbone
    has_3beta_hydroxyl = False
    for match in steroid_matches:
        # Iterate over atoms in the steroid core
        steroid_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
        for atom in steroid_atoms:
            # Check if atom is a carbon atom in a ring
            if atom.GetAtomicNum() == 6 and atom.IsInRing():
                # Check for attached hydroxyl group
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        has_3beta_hydroxyl = True
                        break
            if has_3beta_hydroxyl:
                break
        if has_3beta_hydroxyl:
            break

    if not has_3beta_hydroxyl:
        return False, "No 3β-hydroxyl group found on the steroid backbone"

    return True, "Contains steroid backbone with 3β-hydroxyl group"