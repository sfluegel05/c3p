"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
"""
Classifies: 3beta-hydroxy steroid

Definition: A 3-hydroxy steroid in which the 3-hydroxy substituent is in the beta-position.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # Define SMARTS pattern for steroid backbone
    steroid_smarts = '''
    [#6]1([#6])[#6][#6]2[#6]([#6]1)[#6][#6]3[#6]2[#6][#6][#6]4[#6]3=[#6][#6][#6]([#6]4)
    '''

    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is None:
        return False, "Error in steroid SMARTS pattern"

    # Find matches to the steroid backbone
    matches = mol.GetSubstructMatches(steroid_pattern)
    if not matches:
        return False, "No steroid backbone found"

    # Look for hydroxyl group at position 3
    # For steroids, position 3 is typically the second atom in the pattern
    for match in matches:
        # Get the atom indices from the match
        steroid_atoms = list(match)

        # Assume position 3 is at atom index steroid_atoms[1]
        pos3_idx = steroid_atoms[1]
        pos3_atom = mol.GetAtomWithIdx(pos3_idx)

        # Check if there is a hydroxy (-OH) group attached to position 3
        has_oh = False
        for neighbor in pos3_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                has_oh = True
                break
        if not has_oh:
            continue  # No hydroxyl group at position 3

        # Check if position 3 is a chiral center
        if not pos3_atom.HasProp('_ChiralityPossible'):
            continue  # Position 3 is not chiral

        chiral_tag = pos3_atom.GetChiralTag()
        if chiral_tag == Chem.ChiralType.CHI_UNSPECIFIED:
            continue  # Chirality not specified

        # Get the R/S configuration
        stereo = rdMolDescriptors.CalcAtomRSLabel(pos3_atom)
        if stereo == 'R':
            # In steroids, R configuration at position 3 usually corresponds to beta
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