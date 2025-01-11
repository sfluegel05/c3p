"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion is defined by the presence of a protonated amine group
    and a deprotonated carboxylate group on the same alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino-acid zwitterion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alpha-amino-acid zwitterion pattern: [NH3+][C@H](whatever)[C(=O)[O-]]
    zwitterion_pattern = Chem.MolFromSmarts("[NH3+][C@H]([*!H0])[C](=O)[O-]")
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "Does not match the alpha-amino-acid zwitterion pattern"

    return True, "Matches the alpha-amino-acid zwitterion pattern"