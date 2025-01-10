"""
Classifies: CHEBI:19573 2-enoyl-CoA
"""
from rdkit import Chem

def is_2_enoyl_CoA(smiles: str):
    """
    Determines if a molecule is a 2-enoyl-CoA based on its SMILES string.
    A 2-enoyl-CoA is defined as an unsaturated fatty acyl-CoA in which the S-acyl group contains a double bond between positions 2 and 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-enoyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the specific enoyl-CoA pattern for the double bond at positions 2-3
    enoyl_pattern = Chem.MolFromSmarts("C(=C)C(=O)SCC(N)C(=O)NC(C(=O)O)CCC(=O)NC")  # Capturing a generic enoyl thioester pattern
    
    # Define the CoA moiety pattern using components
    coa_pattern = Chem.MolFromSmarts("COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H](CO[P](=O)(O)O)[C@@H](O)[C@H]1N2C=NC3=C2N=CN=C3N")  # Focusing on the sugar and phosphate backbone with adenine
    
    # Check for the specific enoyl double bond pattern
    if not mol.HasSubstructMatch(enoyl_pattern):
        return False, "No 2-enoyl pattern detected (missing double bond between positions 2 and 3 or thioester linkage)"
    
    # Ensure the Coenzyme A (CoA) structure is present
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A structure detected"

    return True, "Molecule matches the 2-enoyl-CoA definition"