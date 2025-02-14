"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid zwitterion based on its SMILES string.
    An alpha-amino acid zwitterion has a protonated amino group and a deprotonated carboxyl group on the alpha carbon.

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

    # Use a more flexible SMARTS pattern for zwitterionic alpha-amino acid
    zwitterion_pattern = Chem.MolFromSmarts("[C]([NH2+,NH3+])([C](=[O])[O-])")

    matches = mol.GetSubstructMatches(zwitterion_pattern)
    if not matches:
          return False, "Missing zwitterionic alpha-amino acid core structure"

    # Ensure at least one match
    if len(matches) < 1:
        return False, "Missing zwitterionic alpha-amino acid core structure"

    return True, "Molecule contains the zwitterionic alpha-amino acid structure"