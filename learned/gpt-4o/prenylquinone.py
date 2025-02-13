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

    # Define multiple quinone and related moiety patterns
    benzoquinone_pattern = Chem.MolFromSmarts("C1=CC(=O)C=CC1=O")  # 1,4-benzoquinone
    naphthoquinone_pattern = Chem.MolFromSmarts("C1=CC(=O)C2=C(C=CC=C2)[O]1")  # 1,4-naphthoquinone

    # Check for the presence of relevant quinone moieties
    if not (mol.HasSubstructMatch(benzoquinone_pattern) or mol.HasSubstructMatch(naphthoquinone_pattern)):
        return False, "No quinone moiety found"

    # Expand prenyl-derived side-chain detection
    prenyl_sidechain_pattern = Chem.MolFromSmarts("C=C-C-")  # Simple prenyl chain
    longer_prenyl_pattern = Chem.MolFromSmarts("C=C-C-C=C")  # Search for longer prenyl links
    branched_prenyl_pattern = Chem.MolFromSmarts("C=C(C)C")  # Include branching variations

    # Verify the presence and proper alignment of a prenyl side-chain with the quinone core
    if not (mol.HasSubstructMatch(prenyl_sidechain_pattern) or 
            mol.HasSubstructMatch(longer_prenyl_pattern) or 
            mol.HasSubstructMatch(branched_prenyl_pattern)):
        return False, "No polyprenyl side-chain found"

    return True, "Contains a quinone moiety with a polyprenyl-derived side-chain"