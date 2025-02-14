"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone has a flavanone/dihydoflavanone core with a hydroxy substituent at the 3' position
    of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"


    # Combined pattern to identify the flavanone core and 3'-OH
    # [c]1:[c]:[c]:[c](-[C](=O)-[C]2:[c](-[c]:[c]3):[c](:[c]:2):[c]:3)[c]:1 means that we are targeting a carbon in the A-ring with the right substitution pattern, and is connected to a C=O which is then connected to the B-ring
    # [c]1-[c](-[OH])-[c]-[c]-[c]-1 means that we target a 6-membered aromatic ring with an OH in the 2 position.
    # The atom mapping is done so we are sure the 3'-OH group is on the B-ring, using the substructure match.

    flavanone_3_oh_pattern = Chem.MolFromSmarts("[c]1:[c]:[c]:[c](-[C](=O)-[C]2:[c](-[c:1]-[OH]):[c](:[c]:2):[c]:3)[c]:1")


    if not mol.HasSubstructMatch(flavanone_3_oh_pattern):
        return False, "Molecule does not contain a flavanone/dihydoflavanone core with a hydroxyl group at the 3' position"

    return True, "Molecule contains a flavanone/dihydoflavanone core with a hydroxyl group at the 3' position"