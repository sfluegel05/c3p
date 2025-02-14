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

    # 1. Check for carbohydrate moiety (hexose or pentose ring with multiple -OH groups and glycosidic linkages)
    hexose_pattern = Chem.MolFromSmarts("[CX4]([OX2H])[CX4]([OX2H])[CX4]([OX2H])[CX4]([OX2H])[CX4]([OX2H])[CX2H]([OX2H])") # Example hexose
    pentose_pattern = Chem.MolFromSmarts("[CX4]([OX2H])[CX4]([OX2H])[CX4]([OX2H])[CX4]([OX2H])[CX2H2]") # Example pentose
    glycosidic_linkage = Chem.MolFromSmarts("[OX2][CX4]") # C-O-C linkage in carbohydrate

    has_carbohydrate = False
    if mol.HasSubstructMatch(hexose_pattern) or mol.HasSubstructMatch(pentose_pattern):
        carbo_matches = mol.GetSubstructMatches(hexose_pattern)
        carbo_matches.extend(mol.GetSubstructMatches(pentose_pattern))
        has_carbohydrate = True
        if len(mol.GetSubstructMatches(glycosidic_linkage))==0:
          has_carbohydrate=False
    if not has_carbohydrate:
        return False, "No carbohydrate moiety detected (no hexose/pentose ring linked by glycosidic bonds)"

    # 2. Check for lipid moiety (long hydrocarbon chain and lipid-associated functional groups)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([OX1])([OX1])")
    sulfate_pattern = Chem.MolFromSmarts("S(=O)(=O)[OX1]")


    has_lipid = False
    if mol.HasSubstructMatch(long_chain_pattern) or mol.HasSubstructMatch(ester_pattern) or mol.HasSubstructMatch(amide_pattern) or mol.HasSubstructMatch(phosphate_pattern) or mol.HasSubstructMatch(sulfate_pattern):
      has_lipid = True

    if not has_lipid:
        return False, "No lipid moiety detected (no long chain, ester, amide, phosphate or sulfate groups)"
    
    # 3. Check for linkage between carbohydrate and lipid
    if not has_carbohydrate or not has_lipid:
        return False, "No carbohydrate or lipid moiety detected."

    
    #Check number of carbons - must be >10
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Too few carbons for a saccharolipid"

    # 4. Return True and reason
    return True, "Contains both a carbohydrate and a lipid moiety"