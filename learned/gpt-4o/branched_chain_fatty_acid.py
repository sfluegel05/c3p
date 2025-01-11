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
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the carboxylic acid group.
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Define SMARTS pattern for a branch (alkyl substituents like methyl).
    branch_pattern = Chem.MolFromSmarts("C(C)(C)C")
    if mol.HasSubstructMatch(branch_pattern):
        return True, "Contains carboxylic acid group and branched chain (alkyl substituents found)"
    
    # Check for other potential branches, e.g., "-C(C)(C)"
    # To generalize to alkyl branches: Looking for Carbon with 3 single bonds to another carbon
    generic_branches = Chem.MolFromSmarts("C(C)(*)-*")
    matches = mol.HasSubstructMatch(generic_branches)
    
    if matches:
        return True, "Contains carboxylic acid group and branched chain (generic branches found)"
    
    return False, "No suitable alkyl branches found on the hydrocarbon chain"