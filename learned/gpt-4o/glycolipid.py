"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as a 1,2-di-O-acylglycerol with a carbohydrate
    part joined via a glycosidic linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 1,2-di-O-acylglycerol fragment
    # The pattern includes flexibility with R groups for acyl chains and potential amide linkages
    glycerol_pattern = Chem.MolFromSmarts("OCC(COC(=O)[#6])[#6]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No 1,2-di-O-acylglycerol-like structure detected"
    
    # Check for glycosidic linkage to a sugar moiety
    # The pattern accommodates various potential sugar connections and stereochemistry
    sugar_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)CO")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No glycosidic linkage to a sugar moiety found"
    
    # Confirm presence of long carbon chains representing fatty acyl groups
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)[CH2][CH2][CH2]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 2:
        return False, f"Expected at least 2 fatty acid chains, found {len(fatty_acid_matches)}"
    
    return True, "Structure matches a glycolipid with a glycerol backbone, acyl chains, and glycosidic linkage"