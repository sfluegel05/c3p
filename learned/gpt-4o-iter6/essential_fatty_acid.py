"""
Classifies: CHEBI:59549 essential fatty acid
"""
from rdkit import Chem

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is an essential fatty acid based on its SMILES string.
    An essential fatty acid is a polyunsaturated fatty acid required in the diet.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an essential fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for terminal carboxylic acid group (-C(=O)O)
    if not mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)O")):
        return False, "No terminal carboxylic acid group found"
    
    # Count total number of carbons to ensure minimum chain length
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 16:
        return False, f"Insufficient aliphatic chain length (found {carbon_count} carbons), need at least 16"

    # Check for linear chain structure
    linear_chain_pattern = Chem.MolFromSmarts("[#6]-[#6]")
    if not mol.HasSubstructMatch(linear_chain_pattern):
        return False, "Molecule doesn't have a linear hydrocarbon chain"

    # Identify cis double bonds (Z or \), which are critical for polyunsaturation
    cis_double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C/C=C\C|C\C=C/C")))
    if cis_double_bond_count < 2:
        return False, f"Insufficient cis double bonds (found {cis_double_bond_count}), need at least 2 for polyunsaturation"

    return True, "Contains key characteristics of essential fatty acid: carboxylic acid group, multiple cis double bonds, and sufficiently long linear aliphatic chain"