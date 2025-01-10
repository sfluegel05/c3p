"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid consists of a ceramide or sphingoid backbone with a 
    carbohydrate residue attached via a glycosidic linkage.
    
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

    # Patterns for ceramide or sphingoid backbone
    backbone_patterns = [
        Chem.MolFromSmarts("C(=O)N[C@@H](CO)CCCCCCCCCCCCCCCC"),  # Long chain ceramide
        Chem.MolFromSmarts("C[C@@H](N)O[C@H](CO)"),              # Sphingoid base
        Chem.MolFromSmarts("O[C@H](CO)CC(N)")                   # Additional sphingoid variants
    ]
    if not any(mol.HasSubstructMatch(bp) for bp in backbone_patterns):
        return False, "No ceramide or sphingoid backbone detected"

    # Patterns for common carbohydrate moieties
    sugar_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H](O)C1"),  # Hexose ring
        Chem.MolFromSmarts("O[C@H][C@H](O)[C@H](O)CO")                  # Linear hexose
    ]
    if not any(mol.HasSubstructMatch(sp) for sp in sugar_patterns):
        return False, "No carbohydrate moiety detected"

    # Glycosidic linkage pattern
    linkage_patterns = [
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@H]1CO"),  # Common glycosidic link
        Chem.MolFromSmarts("O[C@H]C[C@H](O)C(O)C")                     # Alternative linkage form
    ]
    if not any(mol.HasSubstructMatch(lp) for lp in linkage_patterns):
        return False, "No glycosidic linkage identified"

    # Check for possible sulfation or complex modifications (optional check)
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)O")
    if mol.HasSubstructMatch(sulfate_pattern):
        return True, "Contains sphingoid/ceramide backbone, carbohydrate moiety, glycosidic linkage, and potential sulfation"

    return True, "Contains a sphingoid or ceramide backbone with carbohydrate moiety attached via a glycosidic linkage"