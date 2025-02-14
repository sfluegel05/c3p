"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: 7-hydroxyisoflavones
"""

from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    A 7-hydroxyisoflavone is a hydroxyisoflavone compound having a hydroxy group at the 7-position.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the isoflavone core with atom mapping
    # Atom 7 is labeled with atom map number 7
    isoflavone_core_smarts = """
    [#6]-1:[#6]:[#6]:[#6]:[#6]:[#6]-1
    :[#8]-[#6]-2(=O)-[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-2
    """
    isoflavone_core_smarts = isoflavone_core_smarts.replace('\n', '')
    isoflavone_core = Chem.MolFromSmarts(isoflavone_core_smarts)

    if isoflavone_core is None:
        return False, "Invalid SMARTS pattern for isoflavone core"

    # Find matches of the isoflavone core
    matches = mol.GetSubstructMatches(isoflavone_core)
    if not matches:
        return False, "Molecule does not contain the isoflavone core"

    # Define the SMARTS pattern for hydroxy group at position 7
    # Position 7 corresponds to a specific carbon in the isoflavone core
    hydroxy_smarts = """
    [#6]-1:[#6]:[#6]:[#6]:[#6]:[#6]-1
    :[#8]-[#6]-2(=O)-[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-2
    >>
    [#6]-1:[#6]:[#6]:[#6]:[#6]:[#6]-1
    :[#8]-[#6]-2(=O)-[#6]:[#6]:[#6]:[#6]:[#6]:[#6]-2
    """

    # Label atom positions in the isoflavone core pattern
    isoflavone_core_atoms = isoflavone_core.GetAtoms()
    # Assign atom map numbers to keep track of positions
    for idx, atom in enumerate(isoflavone_core_atoms):
        atom.SetAtomMapNum(idx + 1)

    # Recreate the molecule from SMARTS with atom mapping numbers
    isoflavone_core = Chem.MolFromSmarts(Chem.MolToSmarts(isoflavone_core))

    # Find matches with atom mapping
    matches = mol.GetSubstructMatches(isoflavone_core, useChirality=False)
    if not matches:
        return False, "Molecule does not contain the isoflavone core with correct atom mapping"

    # Position 7 corresponds to atom map number 3 (aromatic carbon next to oxygen)
    position7_atom_map_num = 3  # Adjust based on the SMARTS pattern

    # Check each match for a hydroxy group at position 7
    for match in matches:
        # Get the atom index in the molecule corresponding to position 7
        atom_idx_pos7 = match[position7_atom_map_num - 1]  # Atom map numbers start from 1
        atom_pos7 = mol.GetAtomWithIdx(atom_idx_pos7)

        # Check if there is a hydroxy group attached to position 7
        has_hydroxy = False
        for neighbor in atom_pos7.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Check if oxygen has hydrogen (i.e., is a hydroxy group)
                if neighbor.GetTotalNumHs() > 0:
                    has_hydroxy = True
                    break
        if has_hydroxy:
            return True, "Molecule contains the 7-hydroxyisoflavone core with hydroxy at position 7"

    return False, "Molecule does not have a hydroxy group at position 7 of the isoflavone core"