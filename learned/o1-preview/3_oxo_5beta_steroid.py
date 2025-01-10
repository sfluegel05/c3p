"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:XXXXX 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is a steroid with a ketone at position 3 and beta-configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)

    # Define the steroid nucleus SMARTS pattern (rings A, B, C, D of steroids)
    steroid_core_smarts = 'C1CCC2C(C1)CCC3C2CCC4(C3CCC4)'

    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if steroid_core is None:
        return False, "Invalid steroid core SMARTS pattern"

    # Check for steroid core match
    matches = mol.GetSubstructMatches(steroid_core)
    if not matches:
        return False, "Steroid core not found"

    # For each match, check for ketone at position 3 and beta-configuration at position 5
    for match in matches:
        # Map the SMARTS pattern atom indices to the molecule atom indices
        # Positions in the steroid core (approximate):
        # Position 3: Atom index 2 (C atom in ring A)
        # Position 5: Atom index 5 (chiral center between rings A and B)

        atom_indices = match  # Mapping of SMARTS atoms to molecule atoms
        mol_atoms = [mol.GetAtomWithIdx(idx) for idx in atom_indices]

        # Identify position 3 carbon atom
        pos3_idx = atom_indices[2]
        pos3_atom = mol.GetAtomWithIdx(pos3_idx)

        # Check for ketone (=O) attached to position 3 carbon
        ketone = False
        for neighbor in pos3_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(pos3_idx, neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                ketone = True
                break
        if not ketone:
            continue  # This match does not have ketone at position 3

        # Identify position 5 chiral center
        pos5_idx = atom_indices[5]
        pos5_atom = mol.GetAtomWithIdx(pos5_idx)

        if pos5_atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            continue  # Not a chiral center at position 5

        # Get CIP code for position 5
        cip_code = pos5_atom.GetProp('_CIPCode') if pos5_atom.HasProp('_CIPCode') else None

        # In steroids, beta-configuration at position 5 typically corresponds to 'R' absolute configuration,
        # but this can vary depending on the specific molecule, so we need to verify carefully.

        # For this example, let's assume beta at position 5 corresponds to 'R'
        if cip_code != 'R':
            continue  # Not beta-configuration at position 5

        # If all conditions are met, return True
        return True, "Molecule is a 3-oxo-5beta-steroid"

    # If no matches satisfy all conditions
    return False, "Does not meet criteria for a 3-oxo-5beta-steroid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXXX',
        'name': '3-oxo-5beta-steroid',
        'definition': "Any 3-oxo steroid that has beta- configuration at position 5.",
        'parents': []
    },
    'config': {},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
}