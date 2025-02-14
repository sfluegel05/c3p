"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is a dicarboxylic acid with a hydroxyl group on the alpha carbon
    relative to a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple(bool, str): True and a reason string if the molecule is a 2-hydroxydicarboxylic acid,
                         False and a reason string otherwise, or None, None if input is invalid
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None  # Indicate invalid SMILES

    # 1. Check for two carboxylic acid groups
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if len(carboxyl_matches) != 2:
       return False, f"Molecule does not have 2 carboxyl groups. Found {len(carboxyl_matches)}"

    # 2. Check for a hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) == 0:
        return False, "Molecule does not contain any hydroxyl group."

    #3 Check for Alpha position of Hydroxyl
    alpha_hydroxyl_pattern = Chem.MolFromSmarts("[CX4]([OH])([CX4](=O)O)")
    alpha_hydroxyl_matches = mol.GetSubstructMatches(alpha_hydroxyl_pattern)
    if len(alpha_hydroxyl_matches) == 0:
        return False, "Hydroxyl group is not on the alpha carbon relative to any carboxylic group."

    #4 Success
    return True, "Molecule is a 2-hydroxydicarboxylic acid"