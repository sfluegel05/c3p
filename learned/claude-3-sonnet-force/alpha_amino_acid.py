"""
Classifies: CHEBI:33704 alpha-amino acid
"""
"""
Classifies: CHEBI:33709 alpha-amino acid

An alpha-amino acid is defined as an amino acid in which the amino group is
located on the carbon atom at the position alpha to the carboxy group.
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define alpha-amino acid pattern (allowing additional substituents)
    alpha_amino_acid_pattern = Chem.MolFromSmarts("[N;H1,H2,H3][C;H1,H2,H3]([N;H1,H2,H3,H4])([C;H1,H2,H3])(=[O;H0,H1])")

    # Check if molecule matches the pattern
    matches = mol.GetSubstructMatches(alpha_amino_acid_pattern)

    # Check for zwitterionic forms
    zwitterion_pattern = Chem.MolFromSmarts("[N+;H1,H2,H3,H4][C;H1,H2,H3]([N+;H1,H2,H3,H4])([C;H1,H2,H3])(=[O-1])")
    zwitterion_matches = mol.GetSubstructMatches(zwitterion_pattern)

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_atoms = mol.GetNumAtoms()
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    n_oxygen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    if matches or zwitterion_matches:
        # Check molecular weight (typical range for alpha-amino acids)
        if 60 < mol_wt < 400:
            # Check atom counts (typical ranges for alpha-amino acids)
            if 1 <= n_nitrogen <= 4 and 1 <= n_oxygen <= 6 and 2 <= n_carbon <= 20:
                return True, "Contains an alpha-amino acid substructure and meets additional criteria"
            else:
                return False, "Atom counts outside typical ranges for alpha-amino acids"
        else:
            return False, "Molecular weight outside typical range for alpha-amino acids"

    # Check for cyclic alpha-amino acids
    cyclic_pattern = Chem.MolFromSmarts("[N;H1,H2,H3][C;H1,H2,H3]([N;H1,H2,H3,H4])([C;H1,H2,H3])(=[O;H0,H1])@[ring1,ring2,ring3]")
    cyclic_matches = mol.GetSubstructMatches(cyclic_pattern)

    if cyclic_matches:
        return True, "Contains a cyclic alpha-amino acid substructure"

    return False, "No valid alpha-amino acid substructure found"