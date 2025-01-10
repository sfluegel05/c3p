"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid zwitterion based on its SMILES string.
    An alpha-amino acid zwitterion is characterized by having:
    1. A central (alpha) carbon atom.
    2. A protonated amino group ([NH3+]) connected to this alpha carbon.
    3. A deprotonated carboxylic acid group (C(=O)[O-]) connected to the same alpha carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid zwitterion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for alpha-amino acid zwitterion
    # [C@@H] or [C@H] indicates a chiral alpha carbon, common in amino acids
    # Make sure patterns allow for chirality and potential side chains
    zwitterion_pattern = Chem.MolFromSmarts("[NH3+][C@@H](-*)[C](=O)[O-] | [NH3+][C@H](-*)[C](=O)[O-]")

    # Check if the pattern matches the molecule
    if not mol.HasSubstructMatch(zwitterion_pattern):
        return False, "No alpha-amino acid zwitterion pattern found"
    
    return True, "Molecule is an alpha-amino acid zwitterion with alpha carbon bound to [NH3+] and C(=O)[O-]."

# This function can now be used to determine if a given SMILES belongs to an alpha-amino acid zwitterion class.