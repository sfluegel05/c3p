"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid is a steroid with an alpha-oriented hydroxyl group at the C17 position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule has stereochemistry information
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)

    # General steroid backbone SMARTS with atom mapping
    steroid_smarts = """
    [#6;R1]1-[#6;R1]-[#6;R1]-[#6;R1]-2-[#6;R1]-[#6;R1]-[#6;R1]-3-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-3-[#6;R1]-[#6;R1]-2-[#6;R1]-1
    """

    # Label the C17 atom in the SMARTS pattern
    steroid_smarts = '''
    [C@$](C3)([C@H]4CC[C@@H](O)CC4)C2C1CC[C@H]1C[C@H]2CC3
    '''

    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return False, "Failed to create steroid SMARTS pattern"

    # Perform substructure match with mapping
    matches = mol.GetSubstructMatches(steroid_pattern)
    if not matches:
        return False, "Molecule does not match steroid backbone"

    # Loop over matches to find one with correct stereochemistry
    for match in matches:
        # Map the C17 atom
        atom_indices = list(match)
        c17_idx = atom_indices[16]  # Assuming C17 is at position 16 in the pattern

        c17_atom = mol.GetAtomWithIdx(c17_idx)

        # Check for hydroxyl group at C17
        has_oh = False
        for neighbor in c17_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Check if the oxygen is part of a hydroxyl group
                if len(neighbor.GetNeighbors()) == 1:
                    has_oh = True
                    break

        if not has_oh:
            continue  # Try next match

        # Check stereochemistry at C17 (alpha orientation)
        stereo = c17_atom.GetProp('_CIPCode') if c17_atom.HasProp('_CIPCode') else None
        if stereo == 'S':
            return True, "Molecule is a 17alpha-hydroxy steroid with correct stereochemistry at C17"
        elif stereo == 'R':
            return False, "Hydroxyl group at C17 is beta-oriented"
        else:
            continue  # Try next match

    return False, "No matching stereochemistry found for 17alpha-hydroxy steroid"