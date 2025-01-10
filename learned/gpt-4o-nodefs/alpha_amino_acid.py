"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid has a central alpha carbon atom bonded to an amino group (NH2), a carboxylic acid group (COOH), and a side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for alpha amino acids: [C](N)(C(=O)O)
    # This pattern identifies a carbon bonded to NH, a COO group, and another carbon (side chain R)
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[C;H2,H3](N)(C(=O)[O,H])[!#1]")

    if mol.HasSubstructMatch(alpha_amino_acid_pattern):
        return True, "Structure matches an alpha-amino acid pattern"

    return False, "No match to the alpha-amino acid structure"