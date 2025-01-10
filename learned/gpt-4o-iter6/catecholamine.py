"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine contains a benzene-1,2-diol ring and an aminoethyl side chain.

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

    # SMARTS pattern for catechol structure (benzene-1,2-diol)
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    
    # SMARTS pattern to detect aminoethyl group attached to a benzene core
    aminoethyl_pattern = Chem.MolFromSmarts("NCC")

    # Check for presence of catechol structure
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (benzene-1,2-diol) structure found"

    # Check for the aminoethyl group
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl side chain found"

    # Ensure the correct attachment point
    # Find match indices for catechol and aminoethyl patterns
    catechol_matches = mol.GetSubstructMatches(catechol_pattern)
    aminoethyl_matches = mol.GetSubstructMatches(aminoethyl_pattern)

    # Basic logic to ensure aminoethyl is attached to the catechol 
    # The first carbon of aminoethyl should be one of the carbons in the catechol ring
    for cat_match in catechol_matches:
        for amino_match in aminoethyl_matches:
            if amino_match[1] in cat_match:
                return True, "Contains catechol structure with aminoethyl side chain correctly attached"

    return False, "Catechol and aminoethyl structures do not connect as required"

# Example usage
smiles_example = "CNC[C@H](O)c1ccc(O)c(O)c1"
print(is_catecholamine(smiles_example))