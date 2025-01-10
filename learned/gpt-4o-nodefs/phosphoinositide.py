"""
Classifies: CHEBI:18179 phosphoinositide
"""
from rdkit import Chem

def is_phosphoinositide(smiles: str):
    """
    Determines if a molecule is a phosphoinositide based on its SMILES string.
    A phosphoinositide typically contains a myo-inositol ring with one or more phosphate groups
    and fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphoinositide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for myo-inositol ring
    inositol_pattern = Chem.MolFromSmarts("C1(CO)C(O)C(O)C(O)C(O)O1")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol structure found"

    # Check for phosphate groups (P(=O)(O)(O))
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "Less than one phosphate group found"
    
    # Define pattern for glycerol backbone linking to inositol with ester linkage
    glycerol_with_fatty_acid = Chem.MolFromSmarts("OC[C@H](O)COP(=O)(O)O")
    if not mol.HasSubstructMatch(glycerol_with_fatty_acid):
        return False, "No glycerol backbone with ester-linked phosphate found"

    # Check for long hydrocarbon fatty acid chains (C=12 or more)
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCCCC")  # At least 12 carbon chain
    if len(mol.GetSubstructMatches(long_chain_pattern)) < 2:  # At least two such chains
        return False, "Not enough long hydrocarbon chains found to match fatty acid chains"

    return True, "Contains myo-inositol ring with phosphate group(s) and long fatty acid chains"

# Example usage
# smiles_str = "a SMILES string for a known phosphoinositide"
# is_phosphoinositide(smiles_str)