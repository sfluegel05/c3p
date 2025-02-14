"""
Classifies: CHEBI:192499 anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid pigments, with a benzopyran-4-one core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Core benzopyran-4-one structure (flavone or flavonol) with carbonyl group
    # Using a more general SMARTS that captures flavones and flavonols
    core_pattern = Chem.MolFromSmarts("c1cc(-c2c(O)c(=O)c3ccccc3o2)ccc1 | c1cc(-c2cc(=O)c3ccccc3o2)ccc1 | c1cc(-c2c(=O)c3ccccc3o2)ccc1")
    if not mol.HasSubstructMatch(core_pattern):
         return False, "No benzopyran-4-one core structure found"

    # 2. Check for at least one hydroxyl or methoxy group.
    num_hydroxyls = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OH]")))
    num_methoxys = len(mol.GetSubstructMatches(Chem.MolFromSmarts("OC")))
    if num_hydroxyls == 0 and num_methoxys == 0:
      return False, "No hydroxyl or methoxy groups found"


    return True, "Benzopyran-4-one core structure with hydroxyl and/or methoxy and additional aromatic ring found"