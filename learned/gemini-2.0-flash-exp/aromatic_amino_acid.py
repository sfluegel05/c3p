"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid contains an aromatic ring and an amino acid moiety (amine and carboxylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for aromatic ring (using lowercase 'c' for aromatic carbon)
    aromatic_pattern = Chem.MolFromSmarts("[c]")
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic ring found"
        
    # Check for amino acid moiety (carboxylic acid and amine group, not necessarily alpha)
    # We use a more general pattern that is tolerant to various linkers
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+]~*C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid moiety found"

    # Check if the aromatic ring and amino acid moiety are connected through a linker.
    # The linker can have a few atoms and can be any type. It also allows for a direct bond.
    # This is now more flexible than before.
    combined_pattern1 = Chem.MolFromSmarts("[c]~*~[NX3,NX4+]~*C(=O)O")
    combined_pattern2 = Chem.MolFromSmarts("[c]~[NX3,NX4+]~*C(=O)O") # Amino group directly attached to aromatic ring
    combined_pattern3 = Chem.MolFromSmarts("[c]~*~*~[NX3,NX4+]~*C(=O)O") # Amino separated by more than 1 atom from aromatic
    
    if not (mol.HasSubstructMatch(combined_pattern1) or mol.HasSubstructMatch(combined_pattern2) or mol.HasSubstructMatch(combined_pattern3)):
        return False, "Aromatic ring and amino acid moiety not connected"


    return True, "Contains both an aromatic ring and an amino acid moiety"