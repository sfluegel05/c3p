"""
Classifies: CHEBI:78608 alpha-amino-acid zwitterion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_zwitterion(smiles: str):
    """
    Determines if a molecule is an alpha-amino-acid zwitterion based on its SMILES string.
    An alpha-amino-acid zwitterion has an amino group and a carboxylate group attached to the same
    carbon atom, which is also connected to another carbon atom (alpha position).

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

    # Look for carboxylate group pattern
    carboxylate_pattern = Chem.MolFromSmarts("C([O-])=O")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not carboxylate_matches:
        return False, "No carboxylate group (-COO-) found"

    # Check if amino and carboxylate groups are attached to the same alpha carbon
    for amino_idx in amino_matches:
        amino_atom = mol.GetAtomWithIdx(amino_idx[0])
        for carboxylate_idx in carboxylate_matches:
            carboxylate_carbon = mol.GetAtomWithIdx(carboxylate_idx[0])
            if amino_atom.GetNeighbors()[0].GetIdx() == carboxylate_carbon.GetIdx():
                alpha_carbon = amino_atom.GetNeighbors()[0]
                if len(alpha_carbon.GetNeighbors()) == 3:  # Alpha position
                    # Check chirality if desired
                    # ...
                    return True, "Molecule is an alpha-amino-acid zwitterion"

    return False, "Amino and carboxylate groups are not attached to the same alpha carbon"