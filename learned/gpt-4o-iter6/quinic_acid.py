"""
Classifies: CHEBI:26493 quinic acid
"""
from rdkit import Chem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is quinic acid or its derivative based on its SMILES string.
    
    Quinic acid is a cyclitol carboxylic acid derivative characterized by a cyclohexane ring
    with multiple hydroxyl groups and a carboxylic acid group.
    The molecule may have additional ester-linked modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is quinic acid or its derivative, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # A more flexible SMARTS pattern for quinic acid, allowing variable stereochemistry
    # Basic structure: cyclohexane with at least 3 hydroxyl and one carboxylic acid
    quinic_acid_pattern = Chem.MolFromSmarts("C1([C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C1)C(=O)O")
    
    # Check for quinic acid core
    if not mol.HasSubstructMatch(quinic_acid_pattern):
        return False, "Missing quinic acid backbone structure"
    
    # Additional patterns for ester-linked or acyl-linked modifications
    acyl_pattern = Chem.MolFromSmarts("C(=O)O[C@@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](C1)O")
    if mol.HasSubstructMatch(acyl_pattern):
        return True, "Quinic acid derivative with acyl-linked group(s)"
    
    # Check for any degree of esterification beyond simple esters
    ester_match = Chem.MolFromSmarts("C(=O)O")
    if mol.HasSubstructMatch(ester_match):
        return True, "Quinic acid with esterification"

    # Simple check if any non-hydrogen atom is linked to the possible quinic acid core
    non_h_main_chain = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[!H]C1[C@@H](O)[C@@H](O)[C@@H](O)[C@@H](O)C1C(=O)O")))
    if non_h_main_chain > 0:
        return True, "Quinic acid core with additional non-hydrogen elements"
    
    return True, "Base quinic acid present"