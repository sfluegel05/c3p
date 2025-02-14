"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI:38756 N-acetyl-amino acid
An N-acyl-amino acid that has acetyl as the acyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for amino acid backbone pattern (N-C-C-C=O)
    amino_acid_pattern = Chem.MolFromSmarts("[N;X3][C;X4][C;X4][C;X3](=O)[O;X2]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not amino_acid_matches:
        return False, "No amino acid backbone found"
    
    # Check if the amino nitrogen has an acetyl group attached
    for match in amino_acid_matches:
        amino_N = match[0]
        amino_atom = mol.GetAtomWithIdx(amino_N)
        for neighbor in amino_atom.GetNeighbors():
            if neighbor.GetSmarts() == "CC(=O)":
                return True, "Contains acetyl group attached to amino nitrogen of an amino acid"
    
    return False, "No acetyl group attached to amino nitrogen"