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

    # Look for ester group with a methyl: generic [CH3O][C](=O)
    methyl_ester_pattern = Chem.MolFromSmarts("COC(=O)")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"
    
    # Look for long carbon chains with some flexibility for double bonds or interruptions
    carbon_chain_pattern = Chem.MolFromSmarts("C(CCC)CC")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No long carbon chain found characteristic of fatty acids"
    
    # Ensure that we only have one ester group pattern
    ester_group_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_group_pattern)
    if len(ester_matches) != 1:
        return False, f"Multiple ester groups found, indicating it might not be a simple FAME (found {len(ester_matches)})"
    
    # Re-check for a moderate number of carbon atoms, ensuring it's more than ester-related
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 7:
        return False, f"Too few carbons for a fatty acid methyl ester (found {c_count})"

    return True, "Contains methyl ester group and a suitable carbon chain for a FAME"