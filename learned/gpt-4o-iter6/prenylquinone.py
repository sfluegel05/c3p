"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the quinone core (benzoquinone or naphthoquinone)
    quinone_pattern = Chem.MolFromSmarts("C1=CC(=O)C=C(C=O)C1")
    naphthoquinone_pattern = Chem.MolFromSmarts("C1=CC2=C(C=C1)C(=O)C=CC2=O")

    if not (mol.HasSubstructMatch(quinone_pattern) or mol.HasSubstructMatch(naphthoquinone_pattern)):
        return False, "No quinone structure found"

    # Define SMARTS pattern for polyprenyl chain (repeated isoprene units)
    prenyl_chain_pattern = Chem.MolFromSmarts("(C=C(C)CCC)*")
    
    # Check for polyprenyl side chain
    if not mol.HasSubstructMatch(prenyl_chain_pattern):
        return False, "No prenyl-derived side chain found"

    # Additional classification checks can go here, if necessary

    return True, "Contains quinone structure with prenyl-derived side chain"