"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine typically has a catechol structure with an ethylamine group (2-carbon chain leading to an amine)
    and may include various substitutions.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for broader catechol pattern (benzene with dihydroxy, not restricted to adjacent position)
    catechol_pattern = Chem.MolFromSmarts("c1c(O)cc(O)ccc1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol (dihydroxybenzene) group found"

    # Ethylamine pattern (C-C-N group connected to a catechol),
    # also allowing potential substitutions or variations on ethyl chain
    ethylamine_pattern = Chem.MolFromSmarts("c1ccc(O)c(O)c1CCN")
    if not mol.HasSubstructMatch(ethylamine_pattern):
        ethylamine_pattern_alt = Chem.MolFromSmarts("c1c(O)c(O)ccc1CC[NH2]")
        if not mol.HasSubstructMatch(ethylamine_pattern_alt):
            return False, "No ethylamine group connected to catechol found"

    return True, "Contains catechol group with an ethylamine linkage or variant"