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

    # Define SMARTS patterns for amino group and carboxyl group
    # We expect these to be attached to the same carbon atom
    amino_group = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=O)][$([C;H2,H1;!$(C=O)])]")
    carboxyl_group = Chem.MolFromSmarts("[$(C(=O)[O;H1,-1,-])]")

    # Check if amino and carboxyl groups are attached to the same carbon atom
    matches_a = mol.GetSubstructMatches(amino_group)
    matches_c = mol.GetSubstructMatches(carboxyl_group)
    
    # Look for common carbon indices in both matches (represents α-carbon)
    alpha_carbons = set(m[1] for m in matches_a).intersection(set(m[0] for m in matches_c))

    if not alpha_carbons:
        return False, "Amino and carboxyl groups are not attached to the same α-carbon"

    # If the structure passes all tests, consider it an amino acid
    return True, "Identified as an amino acid with appropriate functional groups"