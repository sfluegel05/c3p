"""
Classifies: CHEBI:26214 porphyrins
"""
"""
Classifies: Porphyrins
"""
from rdkit import Chem

def is_porphyrins(smiles: str):
    """
    Determines if a molecule is a porphyrin based on its SMILES string.
    A porphyrin is characterized by a macrocyclic structure composed of four pyrrole rings connected via methine bridges.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a porphyrin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define porphyrin core SMARTS pattern
    porphyrin_smarts = """
    [nH]1ccc2c1ccc3c2ccc4c3ccc([nH]4)
    """
    # Remove whitespace and newlines from SMARTS
    porphyrin_smarts = porphyrin_smarts.strip().replace('\n', '')
    query_mol = Chem.MolFromSmarts(porphyrin_smarts)
    if query_mol is None:
        return False, "Invalid porphyrin SMARTS pattern"
    
    # Use Substructure Matching
    if mol.HasSubstructMatch(query_mol):
        return True, "Contains porphyrin core macrocyclic structure"
    else:
        return False, "Does not contain porphyrin core macrocyclic structure"