"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid contains both a carbohydrate and a lipid moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for carbohydrate moiety (pyranose/furanose ring with glycosidic linkages)
    pyranose_pattern = Chem.MolFromSmarts("[CX4]1[OX2][CX4][CX4][CX4][CX4]1")
    furanose_pattern = Chem.MolFromSmarts("[CX4]1[OX2][CX4][CX4][CX4]1")
    glycosidic_linkage = Chem.MolFromSmarts("[OX2][CX4]")

    has_carbohydrate = False
    if mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern):
        if len(mol.GetSubstructMatches(glycosidic_linkage)) >= 1:
           has_carbohydrate = True

    if not has_carbohydrate:
      # No carbohydrate detected
        return False, "No carbohydrate moiety detected (no pyranose/furanose ring with glycosidic bond)"


    # 2. Check for lipid moiety (long hydrocarbon chain and lipid-associated functional groups)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]") # at least 7 C's in a row.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([OX1])([OX1])")
    sulfate_pattern = Chem.MolFromSmarts("S(=O)(=O)[OX1]")

    has_lipid = False
    if (mol.HasSubstructMatch(long_chain_pattern) or
        mol.HasSubstructMatch(ester_pattern) or
        mol.HasSubstructMatch(amide_pattern) or
        mol.HasSubstructMatch(phosphate_pattern) or
        mol.HasSubstructMatch(sulfate_pattern)):
         has_lipid = True

    if not has_lipid:
         return False, "No lipid moiety detected (no long chain, ester, amide, phosphate or sulfate groups)"

    # Check number of carbons - must be >10
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Too few carbons for a saccharolipid"

    # Check molecular weight - saccharolipids are typically large
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
       return False, "Molecular weight too low for saccharolipid"

    # 4. Return True and reason
    return True, "Contains both a carbohydrate and a lipid moiety"