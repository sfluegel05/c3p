"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid zwitterion based on its SMILES string.
    
    An alpha-amino acid zwitterion is characterized by having a protonated amino group ([NH3+], [NH2+] or [NH+])
    and a deprotonated carboxylic acid group (C(=O)[O-]) on the same alpha carbon.

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

    # Ensure the structure has an alpha carbon with [NH3+], [NH2+] or [NH+] and C(=O)[O-] groups
    amino_zwitterion_pattern = Chem.MolFromSmarts("[C;R0]([NH3+,[NH2+],[NH+]])[CH,CH2]C(=O)[O-]")
    if not mol.HasSubstructMatch(amino_zwitterion_pattern):
        return False, "No alpha-amino acid zwitterion pattern found"
    
    return True, "Molecule is an alpha-amino acid zwitterion with alpha carbon bound to [NH3+],[NH2+],[NH+] and C(=O)[O-]"

# Example usage:
# smiles = "[NH3+][C@@H](CCC(=O)NCCc1ccc(O)cc1)C([O-])=O"
# print(is_alpha_amino_acid_zwitterion(smiles))