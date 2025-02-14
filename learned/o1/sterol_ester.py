"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:35366 sterol ester
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the steroid nucleus SMARTS pattern (corrected formatting)
    steroid_nucleus_smarts = '[#6;R1]1-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-1-[#6;R1]2-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-2-[#6;R1]3-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-3-[#6;R1]4-[#6;R1]-[#6;R1]-[#6;R1]-[#6;R1]-4'
    steroid_pattern = Chem.MolFromSmarts(steroid_nucleus_smarts)
    if steroid_pattern is None:
        return False, "Invalid steroid SMARTS pattern"

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus found"

    # Template steroid nucleus with atom mapping for positions
    steroid_nucleus_molblock = '''
  RDKit          2D

 34 36  0  0  0  0            999 V2000
   -0.7500    1.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  1
    0.7500    1.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  2
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  3
    0.7500   -1.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  4
   -0.7500   -1.3000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  5
   -1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  6
    2.2500    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  7
    3.7500    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  8
    4.5000    0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0  9
    3.7500   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 10
    2.2500   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 11
    1.5000   -2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 12
    2.2500   -3.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 13
    3.7500   -3.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 14
    4.5000   -2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 15
    3.7500   -0.8000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0 16
    5.7000   -2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0 17
    0.0000   -2.6000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 18
   -2.2500    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 19
   -3.7500    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 20
   -4.5000    0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 21
   -3.7500   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 22
   -2.2500   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 23
   -1.5000   -2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 24
   -2.2500   -3.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 25
   -3.7500   -3.4000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 26
   -4.5000   -2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 27
   -3.7500   -0.8000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0 28
   -5.7000   -2.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0 29
    0.0000    2.6000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0 30
    5.7000    0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 31
    6.4500    2.1000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 32
    6.4500   -0.5000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 33
    7.2000    0.8000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0 34
  1  2  1  0
  2  3  1  0
  3 11  1  0
 11  4  1  0
  4  5  1  0
  5  6  1  0
  6  1  1  0
  3  7  1  0
  7  8  1  0
  8  9  1  0
  9 15  1  0
 15 10  1  0
 10 11  1  0
  9 31  1  0
 31 32  1  0
 32 34  1  0
 34 33  1  0
 33 31  1  0
 10 16  1  0
 15 17  2  0
  5 18  1  0
  6 19  1  0
 19 20  1  0
 20 21  1  0
 21 27  1  0
 27 22  1  0
 22 23  1  0
 23  6  1  0
 22 28  1  0
 27 29  2  0
  2 30  1  0
M  END
'''

    steroid_template = Chem.MolFromMolBlock(steroid_nucleus_molblock, removeHs=False)
    if steroid_template is None:
        return False, "Invalid steroid template"

    # Perform substructure matching with atom mapping
    match = mol.GetSubstructMatch(steroid_template)
    if not match:
        return False, "Steroid nucleus does not match template"

    # Get the atom index corresponding to position 3 (assuming it's atom map '3' in the template)
    atom_map = {}
    for atom in steroid_template.GetAtoms():
        map_num = atom.GetAtomMapNum()
        if map_num:
            atom_map[map_num] = match[atom.GetIdx()]

    position_3 = atom_map.get(3)
    if position_3 is None:
        return False, "Position 3 not found in molecule"

    # Check if the atom at position 3 is connected via an ester linkage
    # Ester pattern: [CX3](=O)[OX2H0][#6]
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2][#6]')
    if ester_pattern is None:
        return False, "Invalid ester SMARTS pattern"

    esters = mol.GetSubstructMatches(ester_pattern)

    # Check if any ester oxygen is connected to position 3
    is_esterified = False
    for ester in esters:
        ester_oxygen = ester[1]
        ester_carbonyl = ester[0]
        if ester_oxygen in [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(position_3).GetNeighbors()]:
            is_esterified = True
            break

    if not is_esterified:
        return False, "Position 3 is not esterified"

    return True, "Sterol ester with esterification at position 3"