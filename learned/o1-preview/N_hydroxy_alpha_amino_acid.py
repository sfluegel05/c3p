"""
Classifies: CHEBI:50760 N-hydroxy-alpha-amino-acid
"""
"""
Classifies: N-hydroxy-alpha-amino-acid
"""

from rdkit import Chem

def is_N_hydroxy_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-hydroxy-alpha-amino-acid based on its SMILES string.
    An N-hydroxy-alpha-amino-acid is an amino acid in which at least one hydrogen attached to the amino group is replaced by a hydroxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-hydroxy-alpha-amino-acid, False otherwise
        str: Reason for classification
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for alpha-amino acid backbone: nitrogen attached to alpha carbon, which is attached to carboxyl group
    amino_acid_pattern = Chem.MolFromSmarts("[N;$([NX3,H2,H1,H0])]C[C](=O)[O]")

    aa_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not aa_matches:
        return False, "No alpha-amino acid backbone found"

    for match in aa_matches:
        nitrogen_idx = match[0]
        nitrogen_atom = mol.GetAtomWithIdx(nitrogen_idx)

        has_n_hydroxy = False

        # Check for N-OH single bond (direct N-hydroxy substitution)
        for neighbor in nitrogen_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                bond = mol.GetBondBetweenAtoms(nitrogen_idx, neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Check if oxygen is hydroxyl (OH)
                    if neighbor.GetTotalDegree() == 1 and neighbor.GetImplicitValence() == 1:
                        has_n_hydroxy = True
                        break

        # Check for N-hydroxyimino group (N=NOH) attached to amino nitrogen
        if not has_n_hydroxy:
            for neighbor in nitrogen_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 7:  # Neighbor nitrogen
                    bond = mol.GetBondBetweenAtoms(nitrogen_idx, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        # Check if neighbor nitrogen has double bond to oxygen (N=O) and single bond to hydroxyl (OH)
                        o_double_bond = False
                        oh_single_bond = False
                        for n2_neighbor in neighbor.GetNeighbors():
                            if n2_neighbor.GetIdx() == nitrogen_idx:
                                continue
                            if n2_neighbor.GetAtomicNum() == 8:  # Oxygen
                                bond_type = mol.GetBondBetweenAtoms(neighbor.GetIdx(), n2_neighbor.GetIdx()).GetBondType()
                                if bond_type == Chem.rdchem.BondType.DOUBLE:
                                    o_double_bond = True
                                elif bond_type == Chem.rdchem.BondType.SINGLE:
                                    if n2_neighbor.GetTotalDegree() == 1 and n2_neighbor.GetImplicitValence() == 1:
                                        oh_single_bond = True
                        if o_double_bond and oh_single_bond:
                            has_n_hydroxy = True
                            break

        if has_n_hydroxy:
            return True, "Contains N-hydroxy-alpha-amino-acid structure"

    return False, "Does not contain N-hydroxy-alpha-amino-acid structure"