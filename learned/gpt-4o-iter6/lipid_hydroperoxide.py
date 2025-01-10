"""
Classifies: CHEBI:61051 lipid hydroperoxide
"""
from rdkit import Chem

def is_lipid_hydroperoxide(smiles: str):
    """
    Determines if a molecule is a lipid hydroperoxide based on its SMILES string.
    A lipid hydroperoxide is any lipid carrying one or more hydroperoxy (-OOH) substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lipid hydroperoxide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for hydroperoxy groups (-OOH)
    hydroperoxy_pattern = Chem.MolFromSmarts("[O][O]")
    hydroperoxy_matches = mol.GetSubstructMatches(hydroperoxy_pattern)
    if not hydroperoxy_matches:
        return False, "No hydroperoxy group found"
    
    # Check for lipid structure - typically long carbon chains, possibly with double bonds
    # We can check for long contiguous carbon chains
    carbon_chain_pattern = Chem.MolFromSmarts("[C]-[C]-[C]-[C]-[C]")  # Example pattern for a chain
    carbon_chain_matches = mol.GetSubstructMatches(carbon_chain_pattern)
    if not carbon_chain_matches:
        return False, "No extended carbon chain typical of lipids found"
    
    # Additional checks for the presence of potential lipid endings like carboxylic acid groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[O]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return True, "Contains hydroperoxy group(s) within a lipid structure"

    # If significant length of carbon chains is observed, it's considered a lipid
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count >= 10:
        return True, "Contains hydroperoxy group(s) with significant carbon chains typical of lipids"

    return False, "Does not fit lipid hydroperoxide profile"