"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine contains a benzene-1,2-diol ring and an aminoethyl side chain or derivatives with such properties.

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

    # Broad SMARTS pattern for catechol structure (benzene with two adjacent hydroxyl groups)
    catechol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)cc1")
    
    # Flexible pattern for aminoethyl group or related side-chain presence on the catechol 
    # and allowing structural variations or substitutions while ensuring intact core structure
    aminoethyl_pattern = Chem.MolFromSmarts("c1(O)c(O)[c|a]c[c|a]c1CCN")  # Allowing for adjacent attachment

    # Check for presence of the catechol structure
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (benzene-1,2-diol) structure found"
    
    # Match for flexible aminoethyl-like attachment to the catechol
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl-like side chain correctly attached to catechol"

    return True, "Contains catechol structure with aminoethyl-like side chain correctly attached"

# Example usage and testing one from the given list
smiles_example = "CNC[C@H](O)c1ccc(O)c(O)c1"  # Example SMILES for adrenaline
print(is_catecholamine(smiles_example))