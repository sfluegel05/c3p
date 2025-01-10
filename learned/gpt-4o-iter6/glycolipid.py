"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined by the presence of a 1,2-di-O-acylglycerol structure 
    linked to a glycosidic saccharide moiety.

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

    # Look for the glycerol backbone pattern (C-C-C with OH groups)
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with hydroxyl groups found"

    # Look for ester group pattern (two ester linkages from the glycerol)
    ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2 for diacyl configuration"

    # Look for glycosidic linkage (oxygen bridge to a carbohydrate)
    glycosidic_pattern = Chem.MolFromSmarts("O[C@H1]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage to carbohydrate detected"

    # Look for sugar moiety (common monosaccharide ring structures)
    sugar_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C(O)C1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No recognizable sugar moiety found"

    return True, "The structure matches the criteria for a glycolipid with glycerol backbone, ester linkages, and glycosidic bonds to sugars"