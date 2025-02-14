"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid is defined as an amino acid where the amino group is located
    on the carbon atom alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"
    except Exception as e:
        return False, f"Error parsing SMILES: {e}"

    # Define SMARTS pattern for alpha-amino acid
    # The pattern matches a carbon atom (alpha carbon) attached to:
    # - A carboxyl group [C(=O)[O-] or C(=O)O]
    # - An amino group [NH2 or NH3+]
    # - Any other substituent [*]
    alpha_amino_acid_pattern = Chem.MolFromSmarts("""
        [NX3+0;H2,H1;!$(N-C=O)]        # Amino nitrogen (not amide), with 1 or 2 hydrogens
        [CX4H]                          # Alpha carbon (tetrahedral carbon with one hydrogen)
        ([#6,#1,#8,#7,#16,#15])         # Any substituent (R group)
        [CX3](=O)[O-,$(O),$(O[H])]      # Carboxy group (acid or carboxylate)
    """)

    if alpha_amino_acid_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check if the molecule matches the alpha-amino acid pattern
    matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)
    if matches:
        return True, "Contains alpha-amino acid structure"
    else:
        return False, "No alpha-amino acid structure found"