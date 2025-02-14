"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is a 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the neoflavonoid core pattern with aryl at position 4.
    # This pattern uses a more specific SMARTS to capture the 1-benzopyran core.
    # The core is composed of fused benzene and a 6 membered ring with oxygen and double bonds.
    # [c] denotes aromatic carbon, [C] denotes a carbon with 4 explicit bonds.
    # The pattern considers 1-benzopyran-2-one (coumarin) as well.
    # [c]1[c](-[CH]=[CH]-O-[CH](-[c]2[c][c][c][c][c]2)-[CH]=[CH]1) also considers the case where we do not have a carbonyl
    # Variations for a dihydro system are also considered.

    core_pattern1 = Chem.MolFromSmarts('c1ccccc1-c2-O-C(-c3ccccc3)=C-C=C2') # 1-benzopyran with phenyl at pos 4
    core_pattern2 = Chem.MolFromSmarts('c1ccccc1-c2-O-C(=O)-C(-c3ccccc3)=C-C2') # 1-benzopyran-2-one (coumarin) with phenyl at pos 4
    core_pattern3 = Chem.MolFromSmarts('c1ccccc1-c2-O-C(-c3ccccc3)-C-C2') # dihydro-1-benzopyran with phenyl at pos 4
    core_pattern4 = Chem.MolFromSmarts('c1ccccc1-c2-O-C(-c3ccccc3)=C-C-C2') # Partially saturated 1-benzopyran
    core_pattern5 = Chem.MolFromSmarts('c1ccccc1-c2-O-C(=O)-C(-c3ccccc3)-C-C2') # Partially saturated coumarin

    if not (mol.HasSubstructMatch(core_pattern1) or mol.HasSubstructMatch(core_pattern2) or mol.HasSubstructMatch(core_pattern3) or mol.HasSubstructMatch(core_pattern4) or mol.HasSubstructMatch(core_pattern5)):
        return False, "Does not contain a 1-benzopyran core with an aryl substituent at position 4."

    return True, "Contains a 1-benzopyran core with an aryl substituent at position 4."