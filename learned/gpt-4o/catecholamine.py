"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine has a 4-(2-Aminoethyl)pyrocatechol core structure with possible substitutions.

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

    # Look for catechol group pattern (benzene ring with 1,2-dihydroxy)
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (1,2-dihydroxybenzene) group found"

    # Look for 2-aminoethyl group attached to the benzene ring
    aminoethyl_pattern = Chem.MolFromSmarts("c1cc([CH2][CH2]N)cc(O)c1O")
    if not mol.HasSubstructMatch(aminoethyl_pattern):
        return False, "No 2-aminoethyl group attached to catechol found"

    return True, "Contains catechol group with a 2-aminoethyl substitution"