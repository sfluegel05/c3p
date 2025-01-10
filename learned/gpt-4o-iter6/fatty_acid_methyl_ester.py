"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester (FAME) based on its SMILES string.
    A FAME is defined as a carboxylic ester obtained by the formal condensation of a fatty acid with methanol.
    
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

    # Look for ester group with a methyl: generic [CH3O][C](=O) verifying connectivity to methanol
    methyl_ester_pattern = Chem.MolFromSmarts("COC(=O)C")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester linkage found"
    
    # Look for a flexible long carbon chain, accounting for more varied structures
    # We use a more generic pattern allowing multiple carbon attachments with potential double bonds
    carbon_chain_pattern = Chem.MolFromSmarts("CCCCCCC")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found characteristic of fatty acids"
    
    # Ensure no multiple ester group patterns, must have exactly one methyl ester
    ester_group_glob_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_group_glob_pattern)
    if len(ester_matches) != 1:
        return False, f"Multiple or no ester groups found (detected {len(ester_matches)} groups)"
    
    # Check the carbon count, ensuring we meet the typical fatty acid methyl ester length
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # Adjusting for larger typical fatty acid chains
        return False, f"Carbon count too low for fatty acid methyl ester (detected {c_count})"
    
    return True, "Contains a distinct methyl ester group and a valid carbon chain for a FAME"