"""
Classifies: CHEBI:33567 catecholamine
"""
from rdkit import Chem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine generally has a catechol (ortho- or para-dihydroxy benzene) structure
    and an aminoalkyl chain with possible variations, such as stereochemistry.

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
    
    # General catechol pattern: ortho- or para-dihydroxy on benzene
    catechol_pattern = Chem.MolFromSmarts("c1cc(O)c(O)cc1 | c1cc(O)ccc(O)c1")
    if not mol.HasSubstructMatch(catechol_pattern):
        return False, "No generalized catechol moiety found"

    # More flexible aminoalkyl chain pattern
    aminoalkyl_pattern = Chem.MolFromSmarts("NCC")
    if not mol.HasSubstructMatch(aminoalkyl_pattern):
        return False, "No general aminoalkyl chain found"

    # Verify there's some common subclass markers or chiral centers
    # Check for stereochemistry in documented common isomers
    if mol.HasSubstructMatch(Chem.MolFromSmarts("N[C@@H](O)")) or mol.HasSubstructMatch(Chem.MolFromSmarts("N[C@H](O)")):
        return True, "Contains generalized catechol structure with flexible aminoalkyl chain and stereochemistry, classified as catecholamine"

    return True, "Contains generalized catechol structure with flexible aminoalkyl chain, classified as catecholamine"