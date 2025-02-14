"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine typically has a catechol structure with a connected alkylamine that can be ethyl or variably substituted.
    The alkylamine group is generally two carbons, ending in an amine, and may include substitutions or stereochemistry.
    
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

    # Look for catechol core pattern (benzene with two adjacent hydroxyl groups)
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No catechol group (benzene-1,2-diol) found"
    
    # Expanded ethylamine-related patterns allowing for substitutions and stereochemistry
    ethylamine_patterns = [
        Chem.MolFromSmarts("c1cc(O)c(O)cc1CCN"),  # Ethylamine direct
        Chem.MolFromSmarts("c1cc(O)c(O)cc1C(C)N"),  # Isoform
        Chem.MolFromSmarts("c1cc(O)c(O)cc1CC[NH2,NH3+]"),  # Protonated forms
        Chem.MolFromSmarts("c1c(O)cc(O)cc1CNC"),  # Alkyl substitution
        Chem.MolFromSmarts("c1cc(O)c(O)cc1C@[CH,H]C[NH2]"),  # Chiral variations
        Chem.MolFromSmarts("c1cc(O)c(O)cc1CX(C)N"),  # General X linked
        # Add more nuanced patterns tailored to false positive characteristics
    ]

    for pattern in ethylamine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains catechol group with an ethylamine linkage or variant"

    return False, "No relevant ethylamine group found connected to catechol"