"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone based on its SMILES string.
    Acrovestone is a polyphenol with an isoflavone core and glycosidic linkages.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acrovestone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for isoflavone core structure (benzopyran-4-one with a phenyl group at the 2-position)
    # More general pattern to account for substitutions
    isoflavone_pattern = Chem.MolFromSmarts("O=C1C(=COc2ccccc2)c3ccccc3O1")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "No isoflavone core structure found"

    # Check for glycosidic linkage (presence of a sugar moiety)
    # More flexible pattern to account for different sugar moieties
    glycosidic_pattern = Chem.MolFromSmarts("[C@H]1O[C@H]([C@H]([C@@H]([C@H]1O)O)O)CO")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage found"

    # Check for polyphenolic nature (multiple hydroxyl groups on aromatic rings)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 3:
        return False, "Not enough hydroxyl groups for a polyphenol"

    # Check molecular weight (acrovestone-like compounds are typically >300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for acrovestone"

    return True, "Contains isoflavone core with glycosidic linkage and polyphenolic features"