"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA has a 3-hydroxy group in the fatty acid chain and is attached to Coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule matches the 3-hydroxy fatty acyl-CoA class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for Coenzyme A end
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)")

    # Define SMARTS pattern for 3-hydroxy chain end, allowing some flexibility
    # considering possibly long hydrophobic fatty acid tails
    hydroxy_pattern = Chem.MolFromSmarts("C[CH](O)C(=O)")

    # Check for both 3-hydroxy group and CoA structure
    if mol.HasSubstructMatch(coa_pattern) and mol.HasSubstructMatch(hydroxy_pattern):
        return True, "Contains the 3-hydroxy fatty acyl-CoA structure"

    return False, "Does not match the 3-hydroxy fatty acyl-CoA pattern"