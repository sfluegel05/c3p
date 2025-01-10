"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define SMARTS patterns for catechol and aminoethyl
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)cc(O)c1")  # benzene ring with two adjacent hydroxyl groups
    aminoethyl_pattern = Chem.MolFromSmarts("CCN")  # represents the aminoethyl side chain (may include variants like CNC etc.)

    # Check for presence of catechol group
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (benzene-1,2-diol) structure found"

    # Check for presence of aminoethyl side chain
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl side chain found"

    return True, "Contains catechol group with aminoethyl side chain"

# Example usage
smiles_example = "CNC[C@H](O)c1ccc(O)c(O)c1"
is_catecholamine(smiles_example)