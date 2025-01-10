"""
Classifies: CHEBI:17297 UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS patterns for components of UDP-sugar
    # Uracil pattern
    uracil_pattern = Chem.MolFromSmarts("n1ccc2[nH]c(=O)n(C)c2=O1")
    
    # Ribose with diphosphate linkage
    diphosphate_ribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](OP(O)(=O)OP(O)(=O)O)C(O)C(O)C1O")

    # Sugar moiety connected through glycosidic bond
    sugar_glycosidic_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O")

    # Check for uracil presence
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "No uracil moiety found"

    # Check for diphosphate linkage with ribose
    if not mol.HasSubstructMatch(diphosphate_ribose_pattern):
        return False, "No diphosphate ribose linkage found"

    # Check for sugar moiety linked via glycosidic bond
    if not mol.HasSubstructMatch(sugar_glycosidic_pattern):
        return False, "No sugar moiety linked via glycosidic bond"
    
    return True, "Contains uracil, diphosphate ribose linkage, and sugar moiety characteristic of UDP-sugars"