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

    # Define essential quinone patterns, considering variations
    quinone_pattern1 = Chem.MolFromSmarts("C1=CC(=O)C=CC1=O")  # 1,4-benzoquinone
    quinone_pattern2 = Chem.MolFromSmarts("C1=CC(=O)C(=C1)OC")  # common in tetrahydrobenzoquinones

    if not (mol.HasSubstructMatch(quinone_pattern1) or mol.HasSubstructMatch(quinone_pattern2)):
        return False, "No quinone moiety found"

    # Extend prenyl side-chain even if more complex and branched (allow variations)
    prenyl_sidechain_pattern = Chem.MolFromSmarts("C=C-C-C=C")  # Prototype for polyprenyl chains
    extended_prenyl_sidechain_pattern = Chem.MolFromSmarts("C=C(C)C")  # Manage branches in prenyl chain

    # Checking for either basic or extended prenyl pattern
    if not (mol.HasSubstructMatch(prenyl_sidechain_pattern) or mol.HasSubstructMatch(extended_prenyl_sidechain_pattern)):
        return False, "No polyprenyl side-chain found"

    return True, "Contains a quinone moiety with a polyprenyl-derived side-chain"