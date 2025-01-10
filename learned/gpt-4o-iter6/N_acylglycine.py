"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine is characterized by an acyl group attached to the nitrogen atom
    of glycine, which typically involves a C(=O)NC structure where the nitrogen is
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

    # Comprehensive pattern for N-acylglycine: C(=O)NC(C)C(=O)O
    # - C(=O): Acyl group attached to nitrogen
    # - N-C-C(=O)O: Amide linkage followed by glycine fragment
    # Patterns capture variability in acyl groups, additional substructures may be considered
    n_acylglycine_pattern = Chem.MolFromSmarts("C(=O)NC(C)C(=O)O")  
    if not mol.HasSubstructMatch(n_acylglycine_pattern):
        return False, "No N-acylglycine structure found in the molecule"
    
    # Additional criteria can be assessed here if necessary
    
    return True, "Molecule contains an N-acylglycine structure"