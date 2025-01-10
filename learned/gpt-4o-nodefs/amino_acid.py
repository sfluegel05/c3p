"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid typically has an amino group and a carboxyl group attached to the same carbon (α-carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for amino group and carboxyl group attached to the α-carbon
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[C;!$(C=O)](N)([*])[C](=O)O")

    if mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return True, "Identified as an amino acid with amino and carboxyl groups attached to the α-carbon"

    return False, "Amino and carboxyl groups are not attached to the same α-carbon"