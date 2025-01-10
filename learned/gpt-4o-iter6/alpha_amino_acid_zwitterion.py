"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid zwitterion based on its SMILES string.
    An alpha-amino acid zwitterion is characterized by having a protonated amino group ([NH3+])
    and a deprotonated carboxylic acid group (C(=O)[O-]) on the same alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid zwitterion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES to get a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for alpha-amino acid zwitterion
    # Note: C(=*) refers to carbonyl oxygen and [O-] carboxylate as connected to same alpha carbon
    zwitterion_pattern = Chem.MolFromSmarts("[C;H1,H2,H3]([NH3+])-[C;H1,H2,H3]-C(=O)[O-]")

    # Check if the pattern matches the molecule
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "No alpha-amino acid zwitterion pattern found"
    
    return True, "Molecule is an alpha-amino acid zwitterion with alpha carbon bound to [NH3+] and C(=O)[O-]"

# Example usage:
# smiles_list = [
#    "[NH3+][C@@H](CCC(=O)NCCc1ccc(O)cc1)C([O-])=O",  # gamma-glutamyltyramine zwitterion
#    "[NH3+][C@@H](CC[SeH])C([O-])=O",  # L-selenohomocysteine zwitterion
#    # Add more test cases here as needed
# ]
# for sm in smiles_list:
#    print(f"{sm}: {is_alpha_amino_acid_zwitterion(sm)}")