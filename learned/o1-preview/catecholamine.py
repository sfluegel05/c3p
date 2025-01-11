"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is defined as 4-(2-aminoethyl)pyrocatechol and derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define catechol moiety pattern (aromatic ring with two adjacent hydroxyl groups)
    catechol_pattern = Chem.MolFromSmarts('c1ccc(O)c(O)c1')
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol moiety found"
    
    # Define aminoethyl chain attached to aromatic ring
    # The pattern matches aromatic carbon connected to a chain of two carbons and a nitrogen
    aminoethyl_pattern = Chem.MolFromSmarts('[c][C][C][N]')
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No aminoethyl chain attached to aromatic ring found"
    
    return True, "Contains catechol moiety with aminoethyl side chain attached to the ring"