"""
Classifies: CHEBI:33856 aromatic amino acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid includes an amino acid structure with an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the amino acid structure: [N][C][C(=O)O]
    amino_acid_pattern = Chem.MolFromSmarts("N[C][C](=O)O")
    if not mol.HasSubstructMatch(amino_acid_pattern):
        return False, "No amino acid functional group found"

    # Look for an aromatic ring (any common aromatic pattern)
    aromatic_ring_pattern = Chem.MolFromSmarts("a")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic ring found"

    # If both patterns are present, it is an aromatic amino acid
    return True, "Molecule contains both an amino acid group and an aromatic ring"

# Example usage:
# is_aromatic_amino_acid("Nc1ccc(cc1)C(=O)O")  # Should return True