"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine contains a catechol (benzene-1,2-diol) ring and an aminoethyl side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extended SMARTS pattern for catechol structure (benzene-1,2-diol)
    catechol_pattern = Chem.MolFromSmarts("c1(c(O)cc(O)c)c1")
    
    # Extended SMARTS pattern for aminoethyl-like structures with potential variability
    # Allows for various aminoethyl variants including substituents
    aminoethyl_pattern = Chem.MolFromSmarts("[NX3][C,c][C,c]")

    # Check for presence of catechol group
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (benzene-1,2-diol) structure found"

    # Check for presence of aminoethyl structural motif
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl-like side chain found"

    # Ensuring correct spatial and logical combination of patterns
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    aminoethyl_matches = mol.GetSubstructMatches(aminoethyl_pattern)

    # Verify combinations
    if not catechol_matches or not aminoethyl_matches:
        return False, "Required structural motifs do not connect correctly"

    return True, "Contains catechol group with aminoethyl-like side chain"

# Example usage
smiles_example = "CNC[C@H](O)c1ccc(O)c(O)c1"
print(is_catecholamine(smiles_example))