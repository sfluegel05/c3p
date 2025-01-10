"""
Classifies: CHEBI:36498 galactosylceramide
"""
from rdkit import Chem

def is_galactosylceramide(smiles: str):
    """
    Determines if a molecule is a galactosylceramide based on its SMILES string.
    A galactosylceramide usually has a galactose moiety linked to a ceramide,
    which includes amide linkages and long fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a galactosylceramide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define improved SMARTS patterns for galactosylceramide
    # Improved pattern to capture beta- and alpha-D-galactosyl variations
    galactose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](O)[C@H](CO)C(O)C1")
    alternate_galactose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](O)[C@@H](CO)C(O)C1")
    
    # This pattern more generically captures the amide linkage in ceramide, including hydroxyamides
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=O)[CX4]")

    # Using a more specific pattern for hydrocarbon chains considering longer ranges
    long_chain_pattern = Chem.MolFromSmarts("C{12,}")
    
    # Match galactose, considering both configurations
    if not (mol.HasSubstructMatch(galactose_pattern) or mol.HasSubstructMatch(alternate_galactose_pattern)):
        return False, "Missing relevant D-galactosyl group"
    
    # Match amide group of ceramide, accounting for hydroxy variants
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "Missing amide linkage characteristic of ceramide"

    # Check for extended hydrocarbon chain features of ceramide
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Missing long-chain characteristic of ceramide"

    return True, "Contains structural motifs of a galactosylceramide: D-galactosyl and ceramide backbone"