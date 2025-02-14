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

    # Flavanone/dihydoflavanone core pattern - less restrictive
    flavanone_core_pattern = Chem.MolFromSmarts("[c]1[c][c][c]2[C]([C]([C]=[O])[C]2)[c]1")
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Molecule does not contain flavanone/dihydoflavanone core"

    # 3'-hydroxy substitution pattern - targets the 3' position with more precision.
    hydroxy_3_pattern = Chem.MolFromSmarts("[c]1-[c](-[OH])-[c](-[c]-[c]1)-[CH]")


    if not mol.HasSubstructMatch(hydroxy_3_pattern):
        return False, "Molecule does not contain a hydroxyl at 3' position"

    return True, "Molecule contains a flavanone/dihydoflavanone core with a hydroxyl group at the 3' position"