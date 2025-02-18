"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion has an amino group and a carboxyl group attached to the same carbon atom,
    with the amino group protonated and the carboxyl group deprotonated.

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

    # Look for amino group pattern
    amino_pattern = Chem.MolFromSmarts("[NH3+]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    if not amino_matches:
        return False, "No amino group (-NH3+) found"

    # Look for carboxyl group pattern
    carboxyl_pattern = Chem.MolFromSmarts("[O-]C=O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group (-COO-) found"

    # Check if amino and carboxyl groups are attached to the same carbon
    for amino_idx in amino_matches:
        amino_atom = mol.GetAtomWithIdx(amino_idx[0])
        for carboxyl_idx in carboxyl_matches:
            carboxyl_carbon = mol.GetAtomWithIdx(carboxyl_idx[0])
            if amino_atom.GetNeighbors()[0].GetIdx() == carboxyl_carbon.GetIdx():
                # Found alpha-amino-acid zwitterion
                return True, "Molecule is an alpha-amino-acid zwitterion"

    return False, "Amino and carboxyl groups are not attached to the same carbon atom"