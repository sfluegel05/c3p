"""
Classifies: CHEBI:24402 glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Broadened ceramide backbone pattern: long chain with amide linkage and OH groups
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H]CO[C@H](O)")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide or sphingoid backbone detected"

    # Patterns to detect more generalized sugar moieties
    # Include broader SMARTS representing a wide range of common saccharides
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H]1(O[C@H](CO)O[C@@H]([C@H]1O)O)"), # Glucose example
        Chem.MolFromSmarts("[C@H]1(O[C@H](O)[C@@H](CO)O[C@H]1O)"), # Galactose example
        Chem.MolFromSmarts("[C@H]1(O[C@H](O)[C@@H](O)[C@H](O)C1)"), # Generic sugar
    ]

    # Check for any sugar pattern match
    sugar_found = any(mol.HasSubstructMatch(sp) for sp in sugar_patterns)
    if not sugar_found:
        return False, "No carbohydrate moiety detected"

    # Identify a glycosidic linkage, indicative of an ether connection with sugars
    glycosidic_linkage_pattern = Chem.MolFromSmarts("[O][C@H]([C@H])O")  # More generic
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "No glycosidic linkage identified"

    return True, "Contains a sphingoid or ceramide backbone with a carbohydrate moiety attached via a glycosidic linkage"