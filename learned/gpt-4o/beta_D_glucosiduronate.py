"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A glucosiduronate is characterized by a beta-D-glucuronic acid moiety deprotonated
    at the carboxyl group typically connected via an O-linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Refined SMARTS for beta-D-glucuronic acid moiety
    # A cyclic six-membered sugar ring with specific stereochemistry
    # Deprotonated carboxylate on C6 position
    beta_d_glucuronic_acid_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H](C1)C(=O)[O-]")
    
    if not mol.HasSubstructMatch(beta_d_glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid moiety found"
    
    # Check for connectivity: ensure there's an O-linkage
    # Pattern for potential O-linkage where the glucuronic acid moiety is linked
    o_linkage_pattern = Chem.MolFromSmarts("[C@H]([C@H](O)[C@@H](O)[C@H](O)C(=O)[O-])O")
    
    if not mol.HasSubstructMatch(o_linkage_pattern):
        return False, "No beta-D-glucuronic O-linkage found"

    # If all checks pass, it is a beta-D-glucosiduronate
    return True, "Contains beta-D-glucuronic acid moiety with proper O-linkage and deprotonation"