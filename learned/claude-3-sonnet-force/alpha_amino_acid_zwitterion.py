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
    nh3_matches = mol.GetSubstructMatches(nh3_pattern)
    if not nh3_matches:
        return False, "No NH3+ group found"

    # Look for COO- group
    coo_pattern = Chem.MolFromSmarts("[O-]C=O")
    coo_matches = mol.GetSubstructMatches(coo_pattern)
    if not coo_matches:
        return False, "No COO- group found"

    # Check connectivity between NH3+ and COO- groups via alpha carbon
    for nh3_idx in nh3_matches:
        nh3_atom = mol.GetAtomWithIdx(nh3_idx)
        for coo_idx in coo_matches:
            coo_atom = mol.GetAtomWithIdx(coo_idx[0])
            bonds = mol.GetBonds()
            for bond in bonds:
                if bond.GetBeginAtomIdx() == nh3_idx and bond.GetEndAtomIdx() == coo_idx[1]:
                    alpha_carbon_idx = bond.GetEndAtomIdx()
                elif bond.GetEndAtomIdx() == nh3_idx and bond.GetBeginAtomIdx() == coo_idx[1]:
                    alpha_carbon_idx = bond.GetBeginAtomIdx()
                else:
                    continue

                # Check for additional substituents on alpha carbon
                alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
                if alpha_carbon.GetTotalNumHs() + len(alpha_carbon.GetNeighbors()) > 3:
                    continue  # Skip if more than two substituents on alpha carbon

                # Check for ring structures
                if mol.GetAtomWithIdx(alpha_carbon_idx).IsInRing():
                    continue  # Skip ring structures for now (can be added later)

                # Check chirality (optional)
                # ...

                return True, "Molecule contains an alpha-amino-acid zwitterion"

    return False, "No alpha-amino-acid zwitterion pattern found"