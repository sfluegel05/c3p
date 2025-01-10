"""
Classifies: CHEBI:26255 prenylquinone
"""
from rdkit import Chem

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    Prenylquinones typically contain a quinone core structure and prenyl-derived side chains.

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

    # Broad range of quinone core patterns
    quinone_patterns = [
        Chem.MolFromSmarts("O=C1C=CC(=O)C=C1"),  # Para-benzoquinone
        Chem.MolFromSmarts("O=C1C=CC=CC(=O)C1"),  # Naphthoquinone
        Chem.MolFromSmarts("O=C1C=CC=C2C=CC=CC2C1=O"),  # Anthraquinone
        Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2C1=O"),  # Anthraquinone alternative
        Chem.MolFromSmarts("O=Cc1ccccc1=O"),  # Simplified quinone structure
    ]

    if not any(mol.HasSubstructMatch(pat) for pat in quinone_patterns):
        return False, "No suitable quinone core structure found"
    
    # Expanded prenyl-like pattern to cover more variations
    prenyl_patterns = [
        Chem.MolFromSmarts("C=C(C)C"),  # Isoprene unit
        Chem.MolFromSmarts("CC=C"),  # Likely prenyl variations
        Chem.MolFromSmarts("C=C([CH2X4,CX4])C"),  # More flexible prenyl unit with variations
        Chem.MolFromSmarts("C=C(C)CO"),  # Hydroxylated prenyl unit
        Chem.MolFromSmarts("C=C(C)OC"),  # Methoxylated prenyl unit
    ]
    
    if not any(mol.HasSubstructMatch(pat) for pat in prenyl_patterns):
        return False, "No prenyl-like side chains found"

    return True, "Contains a broader defined quinone core with prenyl-derived side chain(s)"