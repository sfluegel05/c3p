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

    # Iterate over atoms to find nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:  # Nitrogen
            continue

        nitrogen = atom
        # Check if nitrogen is connected to an alpha carbon
        for neighbor in nitrogen.GetNeighbors():
            if neighbor.GetAtomicNum() != 6:  # Carbon
                continue
            alpha_carbon = neighbor

            # Check if alpha carbon is connected to a carboxyl group
            has_carboxyl = False
            oxygens = []
            for ac_neighbor in alpha_carbon.GetNeighbors():
                if ac_neighbor.GetAtomicNum() == 8:  # Oxygen
                    oxygens.append(ac_neighbor)

            if len(oxygens) >= 2:
                # Check bonds to oxygens
                double_bonded = False
                single_bonded = False
                for oxygen in oxygens:
                    bond = mol.GetBondBetweenAtoms(alpha_carbon.GetIdx(), oxygen.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        double_bonded = True
                    elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        single_bonded = True
                if double_bonded and single_bonded:
                    has_carboxyl = True

            if has_carboxyl:
                # Now check substitutions on nitrogen
                n_has_oh = False
                n_has_noh = False

                # Check for direct N-OH substitutions
                for n_neighbor in nitrogen.GetNeighbors():
                    if n_neighbor.GetIdx() == alpha_carbon.GetIdx():
                        continue  # Skip alpha carbon

                    if n_neighbor.GetAtomicNum() == 8:  # Oxygen
                        # Check if oxygen is hydroxyl group
                        if n_neighbor.GetDegree() == 1:
                            n_has_oh = True

                    elif n_neighbor.GetAtomicNum() == 7:  # Another nitrogen (possible N-hydroxyimino group)
                        n2 = n_neighbor
                        # Check for N-hydroxyimino group connected to amino nitrogen
                        has_n2_double_bond = False
                        has_n2_oh = False
                        for n2_neighbor in n2.GetNeighbors():
                            if n2_neighbor.GetIdx() == nitrogen.GetIdx():
                                continue  # Skip amino nitrogen
                            if n2_neighbor.GetAtomicNum() == 8:
                                bond = mol.GetBondBetweenAtoms(n2.GetIdx(), n2_neighbor.GetIdx())
                                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and n2_neighbor.GetDegree() == 1:
                                    has_n2_oh = True
                            elif n2_neighbor.GetAtomicNum() == 7:
                                bond = mol.GetBondBetweenAtoms(n2.GetIdx(), n2_neighbor.GetIdx())
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    has_n2_double_bond = True
                        if has_n2_double_bond and has_n2_oh:
                            n_has_noh = True

                if n_has_oh or n_has_noh:
                    return True, "Contains N-hydroxy-alpha-amino-acid structure"

    return False, "Does not contain N-hydroxy-alpha-amino-acid structure"