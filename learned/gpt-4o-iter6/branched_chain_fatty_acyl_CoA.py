"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA has a branched-chain fatty acid attached to Coenzyme A via a thioester bond.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define Coenzyme A pattern (simplified for this purpose)
    coa_smarts = "C(=O)SCCNC(=O)CCNC(=O)C(C)(C)C"
    coa_pattern = Chem.MolFromSmarts(coa_smarts)
    
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Define pattern for branching (at least one tertiary carbon, i.e., a carbon bonded to three other carbons)
    branched_pattern = Chem.MolFromSmarts("[CX4]([C])[C][C]")
    
    if not mol.HasSubstructMatch(branched_pattern):
        return False, "No branched-chain in the fatty acid"

    return True, "Contains CoA moiety and branched-chain in the fatty acid"