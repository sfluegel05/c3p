"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is characterized by an ester linkage derived from methanol 
    and a long aliphatic carbon chain representative of fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for methyl ester group
    ester_methyl_pattern = Chem.MolFromSmarts("OC(=O)[C]")
    if not mol.HasSubstructMatch(ester_methyl_pattern):
        return False, "Ester group with methanol (methyl ester) missing"

    # Look for at least one long carbon chain indicative of a fatty acid (e.g., minimum of 6 consecutive carbons)
    long_chain_pattern = Chem.MolFromSmarts("C(CCCCCC)C")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "Missing sufficient long carbon chain for fatty acid"

    # Allow up to two methyl ester groups if the molecule remains a typical fatty acid methyl structure
    ester_matches = mol.GetSubstructMatches(ester_methyl_pattern)
    if not (1 <= len(ester_matches) <= 2):
        return False, f"Found {len(ester_matches)} methyl ester groups, acceptable is 1 or 2"

    return True, "Structure matches fatty acid methyl ester requirements"