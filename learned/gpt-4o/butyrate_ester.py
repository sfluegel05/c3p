"""
Classifies: CHEBI:50477 butyrate ester
"""
from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is defined as any carboxylic ester where the carboxylic acid component is butyric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define improved butyrate ester pattern
    butyrate_patterns = [
        Chem.MolFromSmarts("C1CCC(C(=O)O[*])C1"),  # Linear butanoate ester
        Chem.MolFromSmarts("CC(C)C(=O)O[*]"),  # Recognize branched butyrate
        Chem.MolFromSmarts("C[2H]CC(=O)O[*]"),  # Include isotopic modifications
        Chem.MolFromSmarts("C[13C]CC(=O)O[*]"),  # Include carbon isotopic labeling
    ]

    matches = any(mol.HasSubstructMatch(p) for p in butyrate_patterns)
    if not matches:
        return False, "No butyrate ester group found"

    return True, "Contains a butyrate ester group"