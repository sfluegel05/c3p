"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is a fatty acyl-CoA with a hydroxyl group attached to the third carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the 3-hydroxy group on a fatty acyl chain (R-CO-CH(OH)-R)
    # The pattern represents a carbon chain with a hydroxyl group on the 3rd carbon,
    # and a coenzyme A thioester group (C(=O)SCCNC(=O)...) attached
    pattern = Chem.MolFromSmarts("CC[C@@H](O)C(=O)SCCNC(=O)")
    if mol.HasSubstructMatch(pattern):
        return True, "Contains 3-hydroxy fatty acyl-CoA structure"

    return False, "Does not match the 3-hydroxy fatty acyl-CoA pattern"