"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid
"""

from rdkit import Chem

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl-L-alpha-amino acid based on its SMILES string.
    An N-acyl-L-alpha-amino acid is any L-alpha-amino acid carrying an N-acyl substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acyl-L-alpha-amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for N-acyl-L-alpha-amino acid backbone
    # [C@H]: chiral alpha carbon with L-configuration
    # (NC(=O)R): N-acylated nitrogen
    # C(=O)O: carboxyl group
    pattern = Chem.MolFromSmarts("[C@H](NC(=O)[#6])[#6]C(=O)O")

    if not mol.HasSubstructMatch(pattern):
        return False, "Does not contain N-acyl-L-alpha-amino acid backbone"

    # Find all matches of the pattern in the molecule
    matches = mol.GetSubstructMatches(pattern)
    for match in matches:
        # Atom indices in the match correspond to the pattern atoms
        alpha_carbon_idx = match[0]
        nitrogen_idx = match[1]
        side_chain_idx = match[2]
        carboxyl_carbon_idx = match[3]

        # Check chirality of alpha carbon
        alpha_carbon = mol.GetAtomWithIdx(alpha_carbon_idx)
        chiral_tag = alpha_carbon.GetChiralTag()
        if chiral_tag != Chem.CHI_TETRAHEDRAL_CCW:
            # Not in L-configuration
            continue

        # Check that nitrogen is acylated (bonded to carbonyl carbon)
        nitrogen = mol.GetAtomWithIdx(nitrogen_idx)
        is_acylated = False
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetIdx() != alpha_carbon_idx and neighbor.GetAtomicNum() == 6:
                # Check if this carbon is a carbonyl carbon (C=O)
                has_double_bond_to_oxygen = any(
                    bond.GetBondType() == Chem.BondType.DOUBLE and
                    bond.GetOtherAtom(neighbor).GetAtomicNum() == 8
                    for bond in neighbor.GetBonds()
                    if bond.GetOtherAtom(neighbor) != nitrogen
                )
                if has_double_bond_to_oxygen:
                    is_acylated = True
                    break
        if not is_acylated:
            continue  # Nitrogen is not acylated

        # All criteria met
        return True, "Molecule is an N-acyl-L-alpha-amino acid"

    return False, "Alpha carbon is not in L-configuration or nitrogen is not acylated"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'N-acyl-L-alpha-amino acid',
        'definition': 'Any L-alpha-amino acid carrying an N-acyl substituent.'
    }
}