"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is either a 1,2-di-O-acylglycerol linked to a carbohydrate part
    via a glycosidic bond or a sphingolipid connected to a sugar part.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glycolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Identify glycerol-based lipid structure or sphingolipid structure
    glycerol_pattern = Chem.MolFromSmarts("C(CO)(COC=O)(COC=O)")
    sphingolipid_pattern = Chem.MolFromSmarts("NC(CO)C(O)C=O")
    
    has_lipid_structure = (mol.HasSubstructMatch(glycerol_pattern) or 
                           mol.HasSubstructMatch(sphingolipid_pattern))

    if not has_lipid_structure:
        return False, "No 1,2-di-O-acylglycerol or sphingolipid pattern found"
    
    # Identify glycosidic linkage (sugar joined via oxygen)
    glycosidic_linkage_pattern = Chem.MolFromSmarts("C-O[C@H1,C@H2][O-]")
    if mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return True, "Glycolipid structure identified with lipid and glycosidic linkage"
    
    return False, "No glycosidic linkage found, therefore not a glycolipid"