"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid

Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxy group at position 3 in the beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens
    mol = Chem.AddHs(mol)

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)

    # Define SMARTS pattern for steroid backbone without stereochemistry
    # Atom indices are mapped to identify position 3
    steroid_smarts = """
    [#6;R1]1([#6;R2])[#6;R2][#6;R2]2[#6;R2][#6;R2][#6;R2]3
    [#6;R2][#6;R2][#6;R2][#6;R2][#6;R2][#6;R2]3
    [#6;R2][#6;R2][#6;R2]2[#6;R2][#6;R1]1
    """

    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return False, "Error in steroid SMARTS pattern"

    # Find matches to the steroid backbone
    matches = mol.GetSubstructMatches(steroid_pattern)
    if not matches:
        return False, "No steroid backbone found"

    # Loop over matches to the steroid backbone
    for match in matches:
        # Map position 3 carbon (using approximate numbering)
        # Assuming atom with highest degree in ring A is position 3
        ring_info = mol.GetRingInfo()
        atom_degrees = [(idx, mol.GetAtomWithIdx(idx).GetDegree()) for idx in match]
        sorted_atoms = sorted(atom_degrees, key=lambda x: x[1], reverse=True)
        pos3_idx = sorted_atoms[0][0]
        pos3_atom = mol.GetAtomWithIdx(pos3_idx)

        # Check if there is a hydroxy (-OH) group attached to position 3
        has_oh = False
        for neighbor in pos3_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                has_oh = True
                oxygen_atom = neighbor
                break
        if not has_oh:
            continue  # No hydroxyl group at position 3

        # Check if position 3 is a chiral center
        if not pos3_atom.HasProp('_CIPCode'):
            continue  # Chirality not assigned at position 3

        stereo = pos3_atom.GetProp('_CIPCode')
        if stereo == 'R':
            # In steroids, R configuration at position 3 corresponds to beta
            return True, "3beta-hydroxy steroid identified"
        else:
            return False, "Hydroxy group at position 3 is not in beta configuration"

    return False, "No matching 3beta-hydroxy steroid found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:000000',
        'name': '3beta-hydroxy steroid',
        'definition': 'A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.',
        'parents': [],
    },
    'success': True,
    'message': None
}