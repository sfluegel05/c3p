"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as a lipophilic moiety typically a diacylglycerol or other lipid
    linked to a saccharide part via a glycosidic linkage.

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
    
    # Update: Identify 1,2-di-O-acylglycerol or equivalent lipid structures
    # Improved pattern to also account for various connections of acyl groups
    glycerol_like_pattern = Chem.MolFromSmarts("OCCOC(=O)|OC(*)C(COC(=O)[*])[#6]")
    if mol.HasSubstructMatch(glycerol_like_pattern):
        # Check for glycosidic linkage
        # Flexible pattern to accommodate various configurations and types of sugars
        glycosidic_linkage_pattern = Chem.MolFromSmarts("CO[C@H1,C@H2][O-]")
        if mol.HasSubstructMatch(glycosidic_linkage_pattern):
            return True, "Identified as glycolipid with glycerol-like structure and glycosidic linkage"
    else:
        # Handle glycerol-free glycolipids, often seen in sphingolipids
        sphingolipid_pattern = Chem.MolFromSmarts("NCCO[C@H1,C@H2][O-]")
        if mol.HasSubstructMatch(sphingolipid_pattern):
            return True, "Identified as glycolipid with sphingolipid-like structure"
    
    # If neither pattern has matched, it isn't described well by our current method
    return False, "No discernible glycolipid-like patterns detected"