"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
"""

from rdkit import Chem
from rdkit.Chem import rdqueries

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is a steroid with a ketone (=O) functional group at position 3 of the steroid nucleus.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for the steroid nucleus (cyclopentanoperhydrophenanthrene)
    steroid_nucleus_smarts = """
    [#6]1[#6][#6][#6]2[#6]([#6]1)[#6][#6][#6]3[#6]2[#6][#6][#6]4[#6]3[#6][#6][#6][#6]4
    """
    steroid_nucleus = Chem.MolFromSmarts(steroid_nucleus_smarts)
    if steroid_nucleus is None:
        return False, "Error in steroid nucleus SMARTS pattern"

    # Check if the molecule has the steroid nucleus
    if not mol.HasSubstructMatch(steroid_nucleus):
        return False, "Molecule does not contain steroid nucleus"

    # Find the match of the steroid nucleus to map atom indices
    matches = mol.GetSubstructMatches(steroid_nucleus)
    if not matches:
        return False, "Steroid nucleus not properly matched"

    # Assuming the first match is representative
    match = matches[0]

    # Map the atoms in the steroid nucleus to their positions
    steroid_atoms = { "C1": match[0], "C2": match[1], "C3": match[2], "C4": match[3],
                      "C5": match[4], "C6": match[5], "C7": match[6], "C8": match[7],
                      "C9": match[8], "C10": match[9], "C11": match[10], "C12": match[11],
                      "C13": match[12], "C14": match[13], "C15": match[14], "C16": match[15],
                      "C17": match[16] }

    # Check for a ketone (=O) at position 3
    position_3 = steroid_atoms["C3"]
    atom_3 = mol.GetAtomWithIdx(position_3)
    ketone_found = False
    for neighbor in atom_3.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(position_3, neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
            ketone_found = True
            break

    if not ketone_found:
        return False, "No ketone (=O) group at position 3 of steroid nucleus"

    return True, "Molecule is a 3-oxo steroid"