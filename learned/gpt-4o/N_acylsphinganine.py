"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    An N-acylsphinganine is defined as a ceramide consisting of sphinganine
    in which one of the amino hydrogens is substituted by a fatty acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the sphinganine backbone pattern with flexibility in stereochemistry
    sphinganine_pattern = Chem.MolFromSmarts("[NH][C@H](CO)[C@H](O)CCCCCCCCCCC")
    if not mol.HasSubstructMatch(sphinganine_pattern):
        return False, "No sphinganine backbone found"
    
    # Define the N-acyl pattern with the carbonyl attached to nitrogen
    n_acyl_pattern = Chem.MolFromSmarts("N[C@H](CO)C(=O)C")
    if not mol.HasSubstructMatch(n_acyl_pattern):
        return False, "No N-acyl substitution found"

    # Extend the pattern to include hydroxylation and chain termination in fatty acids
    long_chain_patterns = [
        Chem.MolFromSmarts("C(=O)CCCCCCCCCCC(O)"),  # Hydroxylated fatty acid
        Chem.MolFromSmarts("C(=O)CCCCCCCCCCC"),    # Common fatty acid
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in long_chain_patterns):
        return False, "No suitable long fatty acyl chain connected"

    return True, "Contains sphinganine backbone with N-acyl substitution"