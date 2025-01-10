"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid features a sphingoid or ceramide backbone with a carbohydrate
    residue attached via a glycosidic linkage.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Updated patterns for ceramide or sphingoid backbone
    backbone_patterns = [
        Chem.MolFromSmarts("C(=O)N[C@@H](CO)CCCCCCCCCCCCCCCC"),  # Long chain ceramide
        Chem.MolFromSmarts("C[C@@H](N)O[C@H](CO)"),              # Sphingoid base
        Chem.MolFromSmarts("OC[C@H](NC)C(O)CO")                 # Sphingoid variants
    ]
    if not any(mol.HasSubstructMatch(bp) for bp in backbone_patterns):
        return False, "No ceramide or sphingoid backbone detected"

    # Updated patterns for carbohydrate moieties
    sugar_patterns = [
        Chem.MolFromSmarts("O[C@@H]1[C@@H](O)[C@H](O)[C@H](O)C(O)CO1"),  # Hexose
        Chem.MolFromSmarts("O[C@H]1[C@H](O)C(O)C(O)C(O)C1O"),            # Hexose variant
        Chem.MolFromSmarts("OC[C@@H](O)C1OC(O)C(O)C1O")                 # Pentose
    ]
    if not any(mol.HasSubstructMatch(sp) for sp in sugar_patterns):
        return False, "No carbohydrate moiety detected"

    # Glycosidic linkage detection widened
    linkage_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@@H]1")
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No glycosidic linkage identified"

    # Check for potential sulfation
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")
    if mol.HasSubstructMatch(sulfate_pattern):
        return True, "Contains sphingoid/ceramide backbone, carbohydrate moiety, glycosidic linkage, and sulfate"

    return True, "Contains a sphingoid or ceramide backbone with a carbohydrate moiety attached via a glycosidic linkage"