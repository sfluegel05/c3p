"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is characterized by an acyl group attached to the nitrogen atom
    of glycine, which typically involves a C(=O)NC(C)C(=O)O structure where the nitrogen is
    attached to a -CH2-COOH group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General pattern for N-acylglycine: Look for an amide (N-C(=O)) followed by a glycine-like carboxylate (C-C(=O)O)
    n_acylglycine_pattern = Chem.MolFromSmarts("NC(=O)CC(=O)O")  
    if mol.HasSubstructMatch(n_acylglycine_pattern):
        return True, "Molecule contains an N-acylglycine structure"
    
    # Additional checks or variations could be implemented here as needed
    
    return False, "No N-acylglycine structure found in the molecule"