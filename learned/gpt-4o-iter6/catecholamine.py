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
    
    # SMARTS pattern for aminoethyl group attached to a catechol ring
    # Here, the amine group connects directly to the benzene's side chain at adjacent carbons forming part of the aminoethyl chain
    aminoethyl_pattern = Chem.MolFromSmarts("c1(O)c(O)ccc1CCN")

    # Check for presence of the catechol structure
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (benzene-1,2-diol) structure found"
    
    # Match for aminoethyl attachment to the catechol properly
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl side chain correctly attached to catechol"

    return True, "Contains catechol structure with aminoethyl side chain correctly attached"

# Example usage
smiles_example = "CNC[C@H](O)c1ccc(O)c(O)c1"
print(is_catecholamine(smiles_example))