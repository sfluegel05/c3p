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

    # Look for the glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with hydroxyl groups found"

    # Improved ester group pattern to capture various ester linkages
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Found {len(ester_matches)} ester groups, need at least 2 for diacyl configuration"

    # Enhanced pattern for glycosidic linkages
    glycosidic_pattern = Chem.MolFromSmarts("O[*]")
    if not mol.HasSubstructMatch(glycosidic_pattern):
        return False, "No glycosidic linkage to carbohydrate detected"

    # Broaden pattern for sugar moiety 
    sugar_pattern = Chem.MolFromSmarts("C1OC([CH2,CH,CH]O)C(O)C(O)C1 | C1O[C@H]([C@H](O)[C@H](O)[C@H]1O) | C1OC([CH2,CH,CH]O)C1")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No recognizable sugar moiety found"
    
    return True, "The structure matches the criteria for a glycolipid with glycerol backbone, ester linkages, and glycosidic bonds to sugars"