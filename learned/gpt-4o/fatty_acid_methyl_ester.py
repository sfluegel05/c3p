"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    ester_methyl_pattern = Chem.MolFromSmarts("OC(=O)C")
    if not mol.HasSubstructMatch(ester_methyl_pattern):
        return False, "Ester group with methanol (methyl ester) missing"

    # Count carbon chains; assume fatty acids will have reasonably long carbon chains
    max_long_carbon_chain_length = 0
    for chain in mol.GetSubstructMatches(Chem.MolFromSmarts("[C]!@[C].[C]!@[C]!@[C]!@[C]!@[C]!@[C]")):
        chain_length = len(chain)
        if chain_length > max_long_carbon_chain_length:
            max_long_carbon_chain_length = chain_length

    if max_long_carbon_chain_length < 8:  # strict checks for minimal carbon length for a fatty acid
        return False, "Too few carbon atoms for a fatty acid chain"

    # Ensure exactly one ester linkage derived from methanol is present
    ester_matches = mol.GetSubstructMatches(ester_methyl_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} methyl ester groups, need exactly 1"

    return True, "Structure matches fatty acid methyl ester requirements"