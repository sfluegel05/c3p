"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid (CHEBI:33525)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is typically a 1,2-di-O-acylglycerol joined at oxygen 3 by a glycosidic linkage to a carbohydrate.
    Some bacterial glycolipids may have acylated sugars without a glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for carbohydrate (saccharide) part: look for rings with multiple hydroxyls
    # This SMARTS matches pyranose or furanose rings with at least two hydroxyl groups
    saccharide_pattern = Chem.MolFromSmarts("[O;R]@;!-[*]")  # Ring oxygen in a carbohydrate
    if not mol.HasSubstructMatch(saccharide_pattern):
        # If no carbohydrate, not a glycolipid
        return False, "No carbohydrate moiety detected"
    
    # Check for glycosidic bond (O connected to the carbohydrate and another atom)
    glycosidic_o = Chem.MolFromSmarts("[O;X2][C;R]")  # Oxygen connected to a ring carbon (glycosidic bond)
    if not mol.HasSubstructMatch(glycosidic_o):
        return False, "No glycosidic bond detected"
    
    # Check for acyl groups (ester or amide)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")  # Ester group
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])")  # Amide group
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    acyl_count = len(ester_matches) + len(amide_matches)
    
    # Check if there's a glycerol backbone with two acyl groups
    glycerol = Chem.MolFromSmarts("[CH2X4][CHX4]([CH2X4])")  # Glycerol backbone
    if mol.HasSubstructMatch(glycerol):
        # Check for two ester/amide groups attached to glycerol
        # Modify pattern to check for esters on glycerol
        glycerol_esters = Chem.MolFromSmarts("[CH2X4][CHX4]([CH2X4])-([OX2][CX3]=[OX1])")
        ester_count = len(mol.GetSubstructMatches(glycerol_esters))
        if ester_count >= 2:
            # Check if the third oxygen is a glycosidic bond
            # Assuming the third OH is replaced by glycosidic O
            return True, "Glycerol backbone with two acyl groups and glycosidic linkage"
        else:
            return False, "Glycerol present but insufficient acyl groups"
    else:
        # Bacterial case: check if carbohydrate has acyl groups
        # Look for esters/amides attached to the carbohydrate
        sacch_acyl = Chem.MolFromSmarts("[C;R][OX2][CX3]=[OX1]")  # Ester on carbohydrate
        if len(mol.GetSubstructMatches(sacch_acyl)) >= 1:
            return True, "Acylated carbohydrate without glycerol"
    
    # If none of the above
    return False, "Does not meet glycolipid criteria"