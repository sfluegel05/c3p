"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid is a steroid with a hydroxy group at position 17 in beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the steroid core SMARTS pattern with atom mapping to identify position 17
    steroid_core_smarts = """
    [#6]1[C@]2(CC[C@H]3[C@@H](CC[C@]4([H])C[C@@H](O)[C@]4(CC3)C)[C@@H]2CC1)
    """
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Failed to create steroid core SMARTS pattern"

    # Check if molecule contains the steroid core with correct stereochemistry
    matches = mol.GetSubstructMatches(steroid_core, useChirality=True)
    if not matches:
        return False, "Does not contain steroid core with 17beta-hydroxy group"
    
    # For each match, check if the hydroxy group is at position 17 with beta configuration
    for match in matches:
        atom_indices = match  # Atom indices in the molecule that match the SMARTS pattern
        mol_atom_ids = list(atom_indices)

        # Position 17 corresponds to the atom in the SMARTS pattern where the hydroxy group is attached
        # Assuming atom index 13 in the SMARTS corresponds to position 17
        # Map the SMARTS atom indices to molecule atom indices
        position_17_idx = mol_atom_ids[13]  # Adjust index based on actual pattern

        position_17_atom = mol.GetAtomWithIdx(position_17_idx)
        # Check if the atom at position 17 is a carbon bonded to a hydroxy group
        has_oh = False
        for neighbor in position_17_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                # Check if the oxygen is bonded to hydrogen (hydroxy group)
                for o_neighbor in neighbor.GetNeighbors():
                    if o_neighbor.GetIdx() != position_17_atom.GetIdx() and o_neighbor.GetAtomicNum() == 1:
                        has_oh = True
                        break
        if not has_oh:
            continue  # Try next match

        # Check the stereochemistry at position 17
        if position_17_atom.HasProp('_CIPCode'):
            cip_code = position_17_atom.GetProp('_CIPCode')
            if cip_code == 'R':
                return True, "Found 17beta-hydroxy group with beta configuration"
            else:
                return False, "Hydroxy group at position 17 is not in beta configuration"
        else:
            return False, "Atom at position 17 is not chiral"
    
    # If none of the matches have the 17beta-hydroxy group
    return False, "Did not find 17beta-hydroxy group with beta configuration"