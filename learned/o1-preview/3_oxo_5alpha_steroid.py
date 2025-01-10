"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: CHEBI:XXXX 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid is a steroid with a ketone group at position 3 and alpha configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create a template steroid molecule with atom numbering for positions
    steroid_template = Chem.MolFromMolBlock("""
  Mrv1810 05011906482D          

 17 18  0  0  0  0            999 V2000
    2.8664   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.8664    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1539    0.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4415    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.4415   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.1539   -1.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8664    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.8664   -1.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8664    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7328    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -1.7328   -0.7500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -0.8664   -1.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5992    0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -2.5992   -1.2500    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
   -3.4646    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  6  1  0
  2  3  1  0
  3  4  1  0
  4  5  1  0
  5  6  1  0
  4  7  1  0
  7  8  1  0
  8  9  1  0
  9 10  1  0
  8 11  1  0
 11 12  1  0
 12 13  1  0
 13 14  1  0
 12 15  1  0
 13 16  1  0
 15 17  1  0
 16 17  1  0
M  END
""")
    if steroid_template is None:
        return False, "Error creating steroid template molecule"
    
    # Use the RDKit atom mapping to align the molecule with the template
    pdb_block = Chem.MolToPDBBlock(steroid_template)
    template = Chem.MolFromPDBBlock(pdb_block, removeHs=False)
    if template is None:
        return False, "Error processing steroid template molecule"

    # Generate 3D coordinates for both molecules
    mol = Chem.AddHs(mol)
    template = Chem.AddHs(template)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.EmbedMolecule(template, AllChem.ETKDG())

    # Attempt to align the molecule to the template
    try:
        rmsd = rdMolAlign.GetBestRMS(mol, template)
    except:
        return False, "Could not align molecule with steroid template"

    # Map the molecule to the template to identify positions
    match = mol.GetSubstructMatch(template)
    if not match:
        return False, "Molecule does not match the steroid template"

    # Check for ketone group at position 3
    # Position 3 corresponds to atom index 2 in the template (0-based indexing)
    pos3_atom_idx = match[2]
    pos3_atom = mol.GetAtomWithIdx(pos3_atom_idx)
    # Check if it is a carbonyl carbon (double-bonded to oxygen)
    ketone_found = False
    for neighbor in pos3_atom.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(pos3_atom_idx, neighbor.GetIdx())
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and neighbor.GetAtomicNum() == 8:
            ketone_found = True
            break
    if not ketone_found:
        return False, "No ketone group at position 3"

    # Check for alpha configuration at position 5
    # Position 5 corresponds to atom index 4 in the template
    pos5_atom_idx = match[4]
    pos5_atom = mol.GetAtomWithIdx(pos5_atom_idx)
    # Assign stereochemistry
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    if not pos5_atom.HasProp('_CIPCode'):
        return False, "Chiral center at position 5 is undefined"
    cip_code = pos5_atom.GetProp('_CIPCode')
    if cip_code != 'R':
        return False, "Position 5 does not have alpha configuration"
    
    return True, "Molecule is a 3-oxo-5alpha-steroid"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XXXX',
                              'name': '3-oxo-5alpha-steroid',
                              'definition': 'A 3-oxo steroid that has alpha configuration at position 5.',
                              'parents': []},
        }