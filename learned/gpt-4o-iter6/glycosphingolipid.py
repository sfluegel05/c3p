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

    # Identify the sphingoid base (Ceramide). Typically a long-chain base with an amide.
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](C)CCCCCCCC[C@@H](O)CO") # Example pattern
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone detected"

    # Identify carbohydrate moieties attached via a glycosidic linkage
    sugar_pattern = Chem.MolFromSmarts("[C@H]1(O[C@H]2[C@H]([C@H](O)[C@@H](CO)[C@@H]([C@@H]2O)O)OC1)")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No carbohydrate moiety detected"

    # Check for the glycosidic linkage to O-1 of the sphingoid
    linkage_pattern = Chem.MolFromSmarts("CO[C@H]1") # O-1 linkage pattern
    if not mol.HasSubstructMatch(linkage_pattern):
        return False, "Glycosidic linkage to O-1 not identified"
    
    return True, "Contains a sphingoid base with a carbohydrate moiety attached via glycosidic linkage"