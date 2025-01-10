"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid has an acetyl group attached to the nitrogen of an amino acid backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    found = False
    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen atom
            is_n_acetyl = False
            has_alpha_carbon = False
            has_carboxyl_group = False

            # Check if nitrogen is connected to an acetyl group
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6:  # Carbon atom
                    bonds = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if bonds.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        # Check for carbonyl group (C=O) and methyl group attached to this carbon
                        carbonyl_carbon = nbr
                        double_bonded_oxygen = False
                        methyl_group = False

                        for nbr2 in carbonyl_carbon.GetNeighbors():
                            if nbr2.GetIdx() == atom.GetIdx():
                                continue  # Skip back to nitrogen atom
                            bond = mol.GetBondBetweenAtoms(carbonyl_carbon.GetIdx(), nbr2.GetIdx())
                            if nbr2.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                double_bonded_oxygen = True
                            elif nbr2.GetAtomicNum() == 6 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                # Check if this carbon is a methyl group (attached to 3 hydrogens)
                                num_attached_h = sum(1 for h in nbr2.GetNeighbors() if h.GetAtomicNum() == 1)
                                if num_attached_h == 3:
                                    methyl_group = True

                        if double_bonded_oxygen and methyl_group:
                            is_n_acetyl = True
                            break  # Acetyl group found

            if is_n_acetyl:
                # Now check for the alpha carbon (carbon attached to nitrogen and carboxyl group)
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr != carbonyl_carbon:
                        alpha_carbon = nbr
                        # Check if alpha carbon is connected to a carboxyl group
                        carboxyl_group = False
                        for nbr2 in alpha_carbon.GetNeighbors():
                            if nbr2.GetIdx() == atom.GetIdx():
                                continue  # Skip back to nitrogen atom
                            if nbr2.GetAtomicNum() == 6:
                                # Potential carboxyl carbon
                                carboxyl_carbon = nbr2
                                has_c_double_o = False
                                has_c_single_o = False
                                for nbr3 in carboxyl_carbon.GetNeighbors():
                                    bond = mol.GetBondBetweenAtoms(carboxyl_carbon.GetIdx(), nbr3.GetIdx())
                                    if nbr3.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                        has_c_double_o = True
                                    elif nbr3.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                        has_c_single_o = True
                                if has_c_double_o and has_c_single_o:
                                    carboxyl_group = True
                                    break  # Carboxyl group found
                        if carboxyl_group:
                            has_alpha_carbon = True
                            has_carboxyl_group = True
                            break  # Alpha carbon with carboxyl group found

            if is_n_acetyl and has_alpha_carbon and has_carboxyl_group:
                return True, "Molecule is an N-acetyl-amino acid"

    return False, "Molecule is not an N-acetyl-amino acid"