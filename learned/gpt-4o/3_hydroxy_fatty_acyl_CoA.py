"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA is represented as a fatty acyl-CoA with a hydroxyl group
    attached at the third carbon atom of the fatty acid chain, combined with the CoA moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule matches the 3-hydroxy fatty acyl-CoA class, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the 3-hydroxy fatty acyl-CoA
    # Broaden the capturing of 3rd-hydroxy (tolerant to chain stereochemistry)
    co_a_part = "C(=O)SCCNC(=O)"
    hydroxy_fatty_acid = "CC(O)[CH2:CH]"
    
    # Combine patterns to look for the full motif of coenzyme and terminus
    pattern = Chem.MolFromSmarts(hydroxy_fatty_acid + co_a_part)
    
    # Perform structure matching
    if mol.HasSubstructMatch(pattern):
        return True, "Contains 3-hydroxy fatty acyl-CoA structure"

    return False, "Does not match the 3-hydroxy fatty acyl-CoA pattern"