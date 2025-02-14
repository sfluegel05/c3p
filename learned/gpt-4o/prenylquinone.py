"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define quinone moiety patterns more broadly to capture common variations
    benzoquinone_pattern = Chem.MolFromSmarts("c1cc(=O)c(=O)cc1")  # Generalized 1,4-benzoquinone
    naphthoquinone_pattern = Chem.MolFromSmarts("C1=C(C=CC2=C1C(=O)C(=O)C=C2)C")  # Generalized 1,4-naphthoquinone

    # Check for the presence of relevant quinone moieties
    if not (mol.HasSubstructMatch(benzoquinone_pattern) or mol.HasSubstructMatch(naphthoquinone_pattern)):
        return False, "No quinone moiety found"

    # Enhance prenyl side-chain detection to capture diverse expressions
    simple_prenyl_pattern = Chem.MolFromSmarts("C=C-C")  # Traditional prenyl
    extended_prenyl_pattern = Chem.MolFromSmarts("C=C-C-C=C")  # Common longer prenyl structures
    complex_prenyl_pattern = Chem.MolFromSmarts("CC(C)=C")  # Branching example

    # Ensure the presence and well-positioned prenyl side-chain
    if not (mol.HasSubstructMatch(simple_prenyl_pattern) or 
            mol.HasSubstructMatch(extended_prenyl_pattern) or 
            mol.HasSubstructMatch(complex_prenyl_pattern)):
        return False, "No suitable prenyl side-chain found"

    return True, "Contains a quinone moiety with a polyprenyl-derived side-chain"