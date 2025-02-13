"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.

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

    # Look for ester group pattern with methanol
    ester_methyl_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(ester_methyl_pattern):
        return False, "Ester group with methanol (methyl ester) missing"
    
    # Check for long carbon chains
    # A simple way is to count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Typically, fatty acids have at least 6 carbons in their chain, often many more
    # We include a generic rule that 8 or more is indicative of a fatty acid chain
    if c_count < 8:
        return False, "Too few carbon atoms for a fatty acid"

    # Verify only one ester group is present
    ester_matches = mol.GetSubstructMatches(ester_methyl_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    return True, "Structure matches fatty acid methyl ester requirements"