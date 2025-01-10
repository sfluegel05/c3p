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

    # Accurate SMARTS pattern for catechol (benzene-1,2-diol)
    catechol_pattern = Chem.MolFromSmarts("c1ccc(O)c(O)c1")

    # Accurate SMARTS pattern for aminoethyl group adjacent to the catechol structure
    aminoethyl_pattern = Chem.MolFromSmarts("CCN")
    
    # Check for presence of the catechol structure
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (benzene-1,2-diol) structure found"
    
    # Match for aminoethyl attached to the benzene ring
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No attached aminoethyl group found"

    return True, "Contains catechol structure with an attached aminoethyl group"

# Example usage and testing one from the given list
smiles_example = "CNC[C@H](O)c1ccc(O)c(O)c1"  # Example SMILES for adrenaline
print(is_catecholamine(smiles_example))