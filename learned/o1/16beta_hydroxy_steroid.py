"""
Classifies: CHEBI:17354 16beta-hydroxy steroid
"""
"""
Classifies: CHEBI:XXXXXX 16beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_16beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 16beta-hydroxy steroid based on its SMILES string.
    A 16beta-hydroxy steroid is a steroid with a hydroxy group at position 16 having beta configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 16beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid nucleus template with atom mapping
    steroid_template = """
  Mrv1810 07102312142D          

 21 24  0  0  0  0            999 V2000
    2.0000    0.0000    0.0000 C   0  0  0  0  1  0 16  0  0  0  0  0
    2.8660   -1.5000    0.0000 C   0  0  0  0  2  0 15  0  0  0  0  0
    4.3660   -1.5000    0.0000 C   0  0  0  0  3  0 14  0  0  0  0  0
    5.2321    0.0000    0.0000 C   0  0  0  0  4  0 13  0  0  0  0  0
    4.3660    1.5000    0.0000 C   0  0  0  0  5  0 12  0  0  0  0  0
    2.8660    1.5000    0.0000 C   0  0  0  0  6  0 11  0  0  0  0  0
    2.0000    3.0000    0.0000 C   0  0  0  0  7  0 10  0  0  0  0  0
    0.5000    3.0000    0.0000 C   0  0  0  0  8  0  9  0  0  0  0  0
   -0.3660    1.5000    0.0000 C   0  0  0  0 17  0  8  0  0  0  0  0
    0.5000    0.0000    0.0000 C   0  0  0  0 18  0  1  0  0  0  0  0
    4.3660   -3.0000    0.0000 C   0  0  0  0  9  0 17  0  0  0  0  0
    5.8660   -3.0000    0.0000 C   0  0  0  0 10  0 18  0  0  0  0  0
    6.7321   -1.5000    0.0000 C   0  0  0  0 11  0  4  0  0  0  0  0
    5.8660    0.0000    0.0000 C   0  0  0  0 12  0  4  0  0  0  0  0
    5.8660    1.5000    0.0000 C   0  0  0  0 13  0  5  0  0  0  0  0
    4.3660    3.0000    0.0000 C   0  0  0  0 14  0  5  0  0  0  0  0
   -0.3660   -1.5000    0.0000 C   0  0  0  0 15  0  9  0  0  0  0  0
   -1.8660   -1.5000    0.0000 C   0  0  0  0 16  0  9  0  0  0  0  0
   -2.7321    0.0000    0.0000 C   0  0  0  0 17  0 10  0  0  0  0  0
   -1.8660    1.5000    0.0000 C   0  0  0  0 18  0 11  0  0  0  0  0
    2.8660   -3.0000    0.0000 C   0  0  0  0  8  0 16  0  0  0  0  0
  1  2  1  0
  5  6  1  0
  2  3  2  0
  3  4  1  0
  4  5  2  0
  6  1  1  0
  5 15  1  0
 15 13  1  0
 13 12  2  0
 12 11  1  0
 11  2  1  0
  5 14  1  0
 14 16  1  0
 16  7  1  0
  7  6  1  0
  7  8  1  0
  8  9  2  0
  9 10  1  0
 10  1  1  0
  9 20  1  0
 20 19  1  0
 19 18  2  0
 18 17  1  0
 17 10  1  0
M  END
"""

    # Read the steroid template molecule
    from io import StringIO
    template_mol = Chem.MolFromMolBlock(steroid_template)
    if template_mol is None:
        return False, "Failed to load steroid template"

    # Attempt to align the input molecule to the steroid template
    match = mol.GetSubstructMatch(template_mol)
    if not match:
        return False, "Molecule does not match steroid nucleus"

    # Map the atom indices from the template to the input molecule
    atom_map = {template_idx: mol_idx for template_idx, mol_idx in enumerate(match)}

    # Get the atom corresponding to position 16
    # In the template, let's assume that atom index 15 corresponds to position 16
    # (Note: Atom indices start from 0)
    position_16_atom_idx = atom_map.get(15)  # 15 corresponds to atom index 16 in the template
    if position_16_atom_idx is None:
        return False, "Position 16 not found in molecule"

    position_16_atom = mol.GetAtomWithIdx(position_16_atom_idx)

    # Check if there is a hydroxy group attached to position 16
    is_hydroxy = False
    for neighbor in position_16_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:  # Oxygen atom
            # Check if oxygen is connected to a hydrogen (i.e., it's a hydroxyl group)
            for atom in neighbor.GetNeighbors():
                if atom.GetAtomicNum() == 1:  # Hydrogen atom
                    is_hydroxy = True
                    # Now check stereochemistry at position 16
                    break

    if not is_hydroxy:
        return False, "No hydroxy group at position 16"

    # Check stereochemistry at position 16
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    position_16_chirality = None
    for idx, chirality in chiral_centers:
        if idx == position_16_atom_idx:
            position_16_chirality = chirality
            break

    if position_16_chirality is None:
        return False, "No chiral center at position 16"

    if position_16_chirality != 'R':
        return False, f"Hydroxy group at position 16 does not have beta configuration (chirality: {position_16_chirality})"

    return True, "Molecule is a 16beta-hydroxy steroid"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:XXXXXX',
        'name': '16beta-hydroxy steroid',
        'definition': 'A 16-hydroxy steroid in which the hydroxy group at position 16 has a beta-configuration.',
        'parents': []
    },
    'config': {
        # Configuration details can be added here
    }
}