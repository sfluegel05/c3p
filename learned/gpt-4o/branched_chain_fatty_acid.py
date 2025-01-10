"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a branched-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Define SMARTS pattern for the carboxylic acid group.
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Define a more flexible pattern for branched carbon (branched alkyl side chains).
    # Look for any carbon with more than two connections to other carbons.
    branch_patterns = [
        Chem.MolFromSmarts("C(C)(C)C"),  # simple trimethyl groups
        Chem.MolFromSmarts("C(C)(CC)C"), # slightly longer branches
        Chem.MolFromSmarts("C(C)(C)C=C"), # unsaturated branches
        Chem.MolFromSmarts("C[C@@H](C)O"), # chiral carbon branch
        Chem.MolFromSmarts("[C@H](C)(C)C"), # chiral tertiary carbon
    ]
    
    # Check if at least one branching pattern is found.
    for pattern in branch_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains carboxylic acid group and branched chain (specific pattern matches found)"
    
    # General pattern: a carbon with at least 3 connections, and two separate branches
    generic_branch = Chem.MolFromSmarts("[C&!H0;!$(C=*)](C)(C)[!$(*=*)]") # avoid counting carbonyls as branches
    if mol.HasSubstructMatch(generic_branch):
        return True, "Contains carboxylic acid group and generic branched chain"
    
    return False, "No suitable alkyl branches found on the hydrocarbon chain"