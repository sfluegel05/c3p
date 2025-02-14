"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: chlorophyll

Definition: A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings. The rings can have various side chains which usually include a long phytol chain.
"""
from rdkit import Chem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    A chlorophyll is characterized by a magnesium-containing chlorin ring system
    with a fifth ring and various side chains, often including a long phytol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for chlorophyll core
    # This pattern represents the chlorin ring system with magnesium and the fifth ring
    chlorophyll_smarts = """
    [Mg]
    -[n]1ccc2nc(ccc3nc(ccc4nc(ccc([n]1)c2)c([n]3)c4))-c5ccccc5
    """
    # Remove whitespace and newlines from SMARTS
    chlorophyll_smarts = "".join(chlorophyll_smarts.strip().split())

    # Create RDKit molecule from SMARTS
    chlorophyll_pattern = Chem.MolFromSmarts(chlorophyll_smarts)
    if chlorophyll_pattern is None:
        return False, "Invalid SMARTS pattern for chlorophyll"

    # Perform substructure match
    if mol.HasSubstructMatch(chlorophyll_pattern):
        return True, "Molecule matches chlorophyll core structure"
    else:
        return False, "Molecule does not match chlorophyll core structure"

__metadata__ = {
    'chemical_class': {
        'name': 'chlorophyll',
        'definition': 'A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings. The rings can have various side chains which usually include a long phytol chain.'
    }
}