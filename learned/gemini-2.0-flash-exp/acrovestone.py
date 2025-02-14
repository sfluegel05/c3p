"""
Classifies: CHEBI:2440 acrovestone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acrovestone(smiles: str):
    """
    Determines if a molecule is an acrovestone based on its SMILES string.
    Acrovestones are polyphenols, often with an isoflavone or flavone backbone,
    multiple hydroxyl groups, and often glycosylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acrovestone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for a flavone or isoflavone-like core (more general)
    # Core patterns - allow for variations and dihydro forms
    core_pattern1 = Chem.MolFromSmarts("c1cc(cc(c1)O)C(=O)c2ccccc2O")  # flavone/isoflavone core, some flexibility
    core_pattern2 = Chem.MolFromSmarts("c1cc(cc(c1)O)C(=O)Cc2ccccc2O")  # dihydroflavone
    core_pattern3 = Chem.MolFromSmarts("c1cc(cc(c1)O)C=Cc2ccccc2O")  # flavone core with double bond
    core_pattern4 = Chem.MolFromSmarts("c1cc(cc(c1))C(=O)c2ccccc2") # basic flavone/isoflavone with no OH on B-ring, other positions can have OH
    if not any(mol.HasSubstructMatch(pattern) for pattern in [core_pattern1, core_pattern2, core_pattern3, core_pattern4]):
        return False, "No flavone/isoflavone-like core found"

    # 2. Check for glycosylation (optional)
    glycosidic_bond_pattern = Chem.MolFromSmarts("[C;R]-O-[C;R][OX2]") #C-O-C, and terminal O
    sugar_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
    is_glycosylated = len(sugar_matches) > 0

    # 3. Check for multiple hydroxyl groups or other substitutions
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    methoxy_pattern = Chem.MolFromSmarts("OC") #methoxy group
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern))
    methyl_pattern = Chem.MolFromSmarts("C[CX4]") #methyl group, non-aromatic
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    malonyl_pattern = Chem.MolFromSmarts("C(=O)CC(=O)O") #malonyl group
    malonyl_count = len(mol.GetSubstructMatches(malonyl_pattern))
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O") #carboxyl group
    carboxyl_count = len(mol.GetSubstructMatches(carboxyl_pattern))

    #Conditions for classification. Must have core and one of the following:
    #1. at least 2 hydroxyl groups OR
    #2. at least one glycosylation OR
    #3. at least 1 methoxy and/or malonyl or carboxyl substituent OR
    #4. hydroxyl groups and other substitutions (methyl)
    if hydroxyl_count >= 2 or is_glycosylated or methoxy_count>0 or malonyl_count>0 or carboxyl_count>0 or (hydroxyl_count > 0 and methyl_count>0):
        return True, "Matches acrovestone structural features"

    return False, "Does not match acrovestone characteristics"