"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
"""
Classifies: CHEBI:29412 "alpha-amino-acid zwitterion"
An amino acid-zwitterion obtained by transfer of a proton from the carboxy to the amino group of any alpha-amino acid; major species at pH 7.3.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.

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

    # Look for NH3+ group
    nh3_pattern = Chem.MolFromSmarts("[NH3+]")
    nh3_match = mol.GetSubstructMatch(nh3_pattern)
    if not nh3_match:
        return False, "No NH3+ group found"

    # Look for COO- group
    coo_pattern = Chem.MolFromSmarts("[O-]C=O")
    coo_match = mol.GetSubstructMatch(coo_pattern)
    if not coo_match:
        return False, "No COO- group found"

    # Check connectivity between NH3+ and COO- groups via alpha carbon
    nh3_atom = mol.GetAtomWithIdx(nh3_match[0])
    coo_atom = mol.GetAtomWithIdx(coo_match[1])  # Index 1 is the carbon atom
    alpha_carbon_neighbors = [n.GetIdx() for n in coo_atom.GetNeighbors()]
    if nh3_match[0] in alpha_carbon_neighbors:
        alpha_carbon_idx = coo_match[1]
    else:
        return False, "NH3+ and COO- groups not connected via alpha carbon"

    # Check for additional substituents on alpha carbon
    alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
    if alpha_carbon.GetTotalNumHs() + len(alpha_carbon.GetNeighbors()) > 3:
        return False, "More than two substituents on alpha carbon"

    # Check for ring structures (optional, can be removed if not required)
    if alpha_carbon.IsInRing():
        return False, "Ring structures not considered for this class"

    # Check chirality (optional)
    # ...

    return True, "Molecule contains an alpha-amino-acid zwitterion"