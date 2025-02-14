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
    # This pattern allows for a variety of linkers between the amine and the acid and also accounts for protonated amine.
    amino_acid_pattern = Chem.MolFromSmarts("[NX3,NX4+]~*C(=O)O")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid moiety found"

    # Check if aromatic ring and amino acid moiety are on the same molecule
    # If the two patterns were found separately it does not mean that they are necessarily part of the same molecule.
    # For example, a molecule consisting of benzene and glycine would satisfy the conditions above but is not an aromatic amino acid.
    combined_pattern = Chem.MolFromSmarts("[c]~*~[NX3,NX4+]~*C(=O)O")
    if not mol.HasSubstructMatch(combined_pattern):
      return False, "Aromatic ring and amino acid moiety not connected"
      
    return True, "Contains both an aromatic ring and an amino acid moiety"