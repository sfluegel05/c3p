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

    # Broadened pattern for ceramide: long chain base, typically with an amide linkage & an OH group
    # Simplified pattern to be more inclusive
    ceramide_pattern = Chem.MolFromSmarts("C(=O)N[C@@H]CO")
    if not mol.HasSubstructMatch(ceramide_pattern):
        return False, "No ceramide backbone detected"

    # Patterns to detect common sugar moieties
    # Adjusted SMARTS to include broader configurations of cyclic acetal groups
    sugar_patterns = [
        Chem.MolFromSmarts("[C@H]1(O[C@@H]([C@H](O)[C@@H](CO)O)C(O)C1)"),  # Example glucose
        Chem.MolFromSmarts("[C@@H]1(O[C@H]([C@@H](O)CO)C(O)C(O1)O)")  # Example galactose
    ]

    # Check for any sugar pattern match
    sugar_found = any(mol.HasSubstructMatch(sp) for sp in sugar_patterns)
    if not sugar_found:
        return False, "No carbohydrate moiety detected"

    # Glycosidic linkages - bonds typically occur via oxygen atoms
    # Pattern checks for alcohol or ether linkages indicative of glycosidic bonds
    glycosidic_linkage_pattern = Chem.MolFromSmarts("O[C@H](C)CO")  # Simplified pattern
    if not mol.HasSubstructMatch(glycosidic_linkage_pattern):
        return False, "Glycosidic linkage to O-1 not identified"

    return True, "Contains a sphingoid or ceramide backbone with a carbohydrate moiety attached via glycosidic linkage"