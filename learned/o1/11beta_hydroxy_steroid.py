"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid is any steroid that has a hydroxy group at position 11 with beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid skeleton template with atom mapping (atom numbers correspond to standard numbering)
    steroid_molblock = '''
  Mrv2014 07092014422D          

 21 23  0  0  0  0            999 V2000
    1.2650   -0.7300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5299   -0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7949   -0.7300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    3.7949   -2.1900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.5299   -2.9200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2650   -2.1900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -2.9200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2650   -2.1900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2650   -0.7300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    1.4600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2650    2.1900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5299    1.4600    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5299    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7949   -0.7300    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.7949   -2.1900    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5299   -2.9200    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5299   -4.3800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.2650   -5.1100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -4.3800    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.2650   -5.1100    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  6  1  0
  1 10  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  6  7  1  0
  7  8  1  0
  8  9  1  0
  9 10  1  0
 10 11  1  0
 11 12  1  0
 12 13  1  0
 13 14  1  0
 14  9  1  0
 14 15  1  0
 15 16  1  0
 16 17  1  0
 17  8  1  0
 17 18  1  0
 18 19  1  0
 19 20  1  0
 20  7  1  0
M  END
'''

    # Read the steroid template molecule
    steroid_template = Chem.MolFromMolBlock(steroid_molblock, removeHs=False)
    if steroid_template is None:
        return False, "Failed to parse steroid template molblock"

    # Generate coordinates for the input molecule (needed for chirality analysis)
    AllChem.Compute2DCoords(mol)
    AllChem.Compute2DCoords(steroid_template)

    # Try to find the substructure match
    matcher = Chem.rdFMCS.CreateMCS([steroid_template, mol], ringMatchesRingOnly=True, completeRingsOnly=True)
    if matcher.numAtoms == 0:
        return False, "Does not contain steroid nucleus"

    # Get the atom mapping from the matching
    match = mol.GetSubstructMatch(steroid_template)
    if not match:
        return False, "Does not contain steroid nucleus"

    # Create a mapping from template atom indices to molecule atom indices
    atom_map = {template_idx: mol_idx for template_idx, mol_idx in enumerate(match)}

    # Atom numbering in the template corresponds to standard steroid numbering
    # Position 11 corresponds to atom index 10 (since atom indices start from 0)

    atom_11_template_idx = 10  # atom numbering in molfile starts from 1, so index 10 corresponds to atom 11

    # Get the atom index of atom 11 in the molecule
    atom_11_mol_idx = atom_map.get(atom_11_template_idx)
    if atom_11_mol_idx is None:
        return False, "Failed to map atom 11"

    # Get the atom at position 11
    atom_11 = mol.GetAtomWithIdx(atom_11_mol_idx)

    # Check if atom 11 is a carbon
    if atom_11.GetAtomicNum() != 6:
        return False, "Atom at position 11 is not carbon"

    # Check if atom 11 has an attached hydroxyl group (O-H)
    has_oh = False
    for neighbor in atom_11.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:  # Oxygen
            # Check if the oxygen is connected to a hydrogen (i.e., it's a hydroxyl group)
            if neighbor.GetTotalNumHs() > 0:
                has_oh = True
                oxygen_atom = neighbor
                break

    if not has_oh:
        return False, "No hydroxyl group at position 11"

    # Check the configuration at atom 11
    chiral_tag = atom_11.GetChiralTag()
    if chiral_tag == Chem.ChiralType.CHI_UNSPECIFIED:
        return False, "Atom at position 11 is not chiral"

    # Assign stereochemistry (needed for correct assignment)
    Chem.AssignAtomChiralTagsFromStructure(mol, replaceExistingTags=True)

    # Get the CIP code at atom 11
    try:
        cip_code = atom_11.GetProp('_CIPCode')
    except KeyError:
        return False, "Cannot determine stereochemistry at position 11"

    if cip_code != 'R':
        return False, f"Hydroxy group at position 11 is not beta (CIP code: {cip_code})"

    return True, "Contains 11beta-hydroxy group on steroid nucleus"