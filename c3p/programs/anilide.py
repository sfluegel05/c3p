"""
Classifies: CHEBI:13248 anilide
"""
"""
Classifies: anilide
"""
from rdkit import Chem

def is_anilide(smiles: str):
    """
    Determines if a molecule is an anilide based on its SMILES string.
    An anilide is any aromatic amide obtained by acylation of aniline,
    meaning it has an amide group where the nitrogen is directly attached
    to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anilide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for anilide
    # Pattern explanation:
    # - [a]: Any aromatic atom (includes heteroatoms in aromatic rings)
    # - [NX3]: Nitrogen atom with three connections (amine nitrogen)
    # - [CX3](=O): Carbonyl carbon (double-bonded to oxygen)
    anilide_pattern = Chem.MolFromSmarts("[a][NX3][CX3](=O)")
    
    if mol.HasSubstructMatch(anilide_pattern):
        return True, "Contains an anilide moiety (amide nitrogen attached to an aromatic ring)"
    else:
        return False, "Does not contain an anilide moiety"