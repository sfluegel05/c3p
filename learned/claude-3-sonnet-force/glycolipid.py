"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    
    A glycolipid is defined as:
    Any member of a class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage to a carbohydrate part 
    (usually a mono-, di- or tri-saccharide). Some substances classified as bacterial glycolipids have the sugar part 
    acylated by one or more fatty acids and the glycerol part may be absent.

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
    
    # Look for saccharide pattern
    saccharide_pattern = Chem.MolFromSmarts("[OX2][C;!$(C=O)][OX2]")
    if not mol.HasSubstructMatch(saccharide_pattern):
        return False, "No saccharide pattern found"
    
    # Look for glycosidic linkage (-O-C-O-C-)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[OX2][C;!$(C=O)][OX2][C;!$(C=O)]")
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "No glycosidic linkage found"
    
    # Look for fatty acid chains or acylated saccharides
    acyl_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    acyl_matches = mol.GetSubstructMatches(acyl_pattern)
    if not acyl_matches:
        return False, "No fatty acid chains or acylated saccharides found"
    
    # Check for glycerol backbone (optional for bacterial glycolipids)
    glycerol_pattern = Chem.MolFromSmarts("[OX2][C;!$(C=O)][C;!$(C=O)][OX2]")
    if mol.HasSubstructMatch(glycerol_pattern):
        return True, "Contains glycerol backbone with glycosidic linkage and acyl groups"
    else:
        return True, "Bacterial glycolipid with acylated saccharide"
    
    # Additional checks/filters can be added here if needed