"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
         bool: True if molecule is a 4'-hydroxyflavanone, False otherwise.
         str:  Reason for the classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Combined SMARTS pattern for flavanone core with 4'-hydroxy substitution
    # Numbering is as follows: A ring (left) is c1ccc2...; C ring is ..C(=O)CC..; and B ring (right) is ..(Oc2c1)c1ccccc1
    # The 4' carbon is the 4th carbon in this ring, starting from the link to the rest of the molecule.
    # Note that the numbering is not the canonical one in the documentation.
    # Using explicit H atoms on the carbons, otherwise we match anything that is bonded to 4 things
    # and we do not want to match the carbonyl carbon
    combined_pattern = Chem.MolFromSmarts("c1[cH]c[cH]c2[C](=O)[CH2][CH](Oc2c1)[cH]1[cH][cH][c]([OH])[cH][cH]1")

    if mol.HasSubstructMatch(combined_pattern):
         return True, "Molecule contains a flavanone core with a hydroxy substituent at the 4' position of the phenyl ring"
    else:
        return False, "Molecule does not have a hydroxyl substituent at the 4' position of the phenyl ring"