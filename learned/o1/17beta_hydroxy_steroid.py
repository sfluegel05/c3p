"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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
    steroid_core_smarts = '[#6@H]1[C@@H]2CC[C@H]3[C@@H](CC[C@]4([H])[C@H](O)[C@@H](C[C@H]4C3)[C@]2([H])CC1)[C@]([H])(C)C'
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Failed to create steroid core SMARTS pattern"

    # Check if molecule contains the steroid core with correct stereochemistry
    matches = mol.GetSubstructMatches(steroid_core, useChirality=True)
    if not matches:
        return False, "Does not contain steroid core with 17beta-hydroxy group"

    # For each match, check if the hydroxy group is at position 17 with beta configuration
    for match in matches:
        # Map the atom indices
        # The SMARTS pattern is designed so that the atom connected to the OH group is at a specific position
        match_atoms = {i: mol.GetAtomWithIdx(idx) for i, idx in enumerate(match)}
        
        # Identify the atom corresponding to position 17 (e.g., atom index 7 in the SMARTS pattern)
        # Adjust the index based on the actual SMARTS pattern
        position_17_atom = match_atoms[7]

        # Check if the atom at position 17 is a carbon bonded to a hydroxy group
        has_oh = False
        for neighbor in position_17_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen atom
                # Check if the oxygen is bonded to hydrogen (hydroxy group)
                for o_neighbor in neighbor.GetNeighbors():
                    if o_neighbor.GetIdx() != position_17_atom.GetIdx() and o_neighbor.GetAtomicNum() == 1:
                        has_oh = True
                        break
                if has_oh:
                    break
        if not has_oh:
            continue  # Try next match

        # Check the stereochemistry at position 17
        # In RDKit, chiral centers can be assigned R/S configuration
        # We'll use CIP codes to determine the configuration at position 17
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        if position_17_atom.HasProp('_CIPCode'):
            cip_code = position_17_atom.GetProp('_CIPCode')
            if cip_code == 'R':
                # For steroids, position 17 R configuration corresponds to beta orientation
                return True, "Found 17beta-hydroxy group with beta configuration"
            else:
                return False, "Hydroxy group at position 17 is not in beta configuration"
        else:
            return False, "Atom at position 17 is not chiral"
        
    # If none of the matches have the 17beta-hydroxy group
    return False, "Did not find 17beta-hydroxy group with beta configuration"