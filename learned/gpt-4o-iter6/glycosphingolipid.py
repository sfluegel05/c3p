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

    # Enhanced ceramide or sphingoid backbone pattern
    backbone_patterns = [
        Chem.MolFromSmarts("C(=O)N[C@@H](CO)CO"),   # More general backbone
        Chem.MolFromSmarts("C(N[C@H](CO)O)CCCC")    # Include long chain
    ]
    if not any(mol.HasSubstructMatch(bp) for bp in backbone_patterns):
        return False, "No ceramide or sphingoid backbone detected"

    # Broadened carbohydrate detection patterns
    sugar_patterns = [
        Chem.MolFromSmarts("OC[C@H]1O[C@@H](CO)[C@H](O)[C@H](O)[C@H]1O"),  # Generic hexose
        Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"),   # Variation
        # More cyclic ether patterns to generically represent sugars
        Chem.MolFromSmarts("O[C@H]1C[C@H](O)[C@@H](CO)O[C@@H]1O")         # Cyclic ether
    ]
    if not any(mol.HasSubstructMatch(sp) for sp in sugar_patterns):
        return False, "No carbohydrate moiety detected"

    # Enhanced glycosidic linkage pattern
    linkage_pattern = Chem.MolFromSmarts("[C@H]([O])[C@H]([O])")  # Capture common linkage
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "No glycosidic linkage identified"

    return True, "Contains a sphingoid or ceramide backbone with a carbohydrate moiety attached via a glycosidic linkage"