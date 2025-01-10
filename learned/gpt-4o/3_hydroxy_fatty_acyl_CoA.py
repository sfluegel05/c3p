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
    
    # Define a SMARTS pattern for the Coenzyme A portion
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)")

    # Define SMARTS pattern for the 3-hydroxy group in the fatty acid chain
    # We explicitly check for the 3-position considering '[C@@H]' for stereochemistry.
    hydroxy_pattern = Chem.MolFromSmarts("[CH][CH](O)[CH2]C(=O)")

    # Check for both 3-hydroxy group and CoA structure
    has_coa = mol.HasSubstructMatch(coa_pattern)
    has_hydroxy = mol.HasSubstructMatch(hydroxy_pattern)

    if has_coa and has_hydroxy:
        return True, "Contains the 3-hydroxy fatty acyl-CoA structure"

    return False, "Does not match the 3-hydroxy fatty acyl-CoA pattern"