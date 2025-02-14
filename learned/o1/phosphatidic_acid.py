"""
Classifies: CHEBI:16337 phosphatidic acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid is a derivative of glycerol where one hydroxyl group is esterified with phosphoric acid,
    and the other two are esterified with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern with open sites for substitutions
    glycerol_pattern = Chem.MolFromSmarts("[C;H2]-[C;H]-[C;H2]")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "Glycerol backbone not found"

    for match in matches:
        C1_idx, C2_idx, C3_idx = match
        C1 = mol.GetAtomWithIdx(C1_idx)
        C2 = mol.GetAtomWithIdx(C2_idx)
        C3 = mol.GetAtomWithIdx(C3_idx)

        # Check connections for C1 and C3 (esterified fatty acids)
        esterified_fatty_acids = 0
        for idx in [C1_idx, C3_idx]:
            atom = mol.GetAtomWithIdx(idx)
            ester_found = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8:  # Oxygen
                    for nbr in neighbor.GetNeighbors():
                        if nbr.GetIdx() != atom.GetIdx() and nbr.GetAtomicNum() == 6:  # Carbonyl carbon
                            # Check if this carbon is double bonded to oxygen (C=O)
                            carbonyl_carbon = nbr
                            has_double_bond_to_O = False
                            for bond in carbonyl_carbon.GetBonds():
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    other_atom = bond.GetOtherAtom(carbonyl_carbon)
                                    if other_atom.GetAtomicNum() == 8:
                                        has_double_bond_to_O = True
                            if has_double_bond_to_O:
                                ester_found = True
                                esterified_fatty_acids += 1
                                break
                    if ester_found:
                        break
            if not ester_found:
                break  # If one of the carbons doesn't have an ester linkage, move to the next match

        if esterified_fatty_acids != 2:
            continue  # Not all esterified fatty acids found, check next glycerol backbone

        # Check connection for C2 (phosphate group)
        phosphate_found = False
        for neighbor in C2.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for nbr in neighbor.GetNeighbors():
                    if nbr.GetAtomicNum() == 15:  # Phosphorus
                        phosphorus = nbr
                        # Check for P(=O)(O)(O)
                        dbl_bonded_O = 0
                        single_bonded_O = 0
                        for bond in phosphorus.GetBonds():
                            other = bond.GetOtherAtom(phosphorus)
                            if other.GetAtomicNum() == 8:  # Oxygen
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    dbl_bonded_O += 1
                                elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                                    single_bonded_O += 1
                        if dbl_bonded_O >= 1 and single_bonded_O >= 2:
                            phosphate_found = True
                            break
                if phosphate_found:
                    break

        if not phosphate_found:
            continue  # Phosphate group not found at C2, check next glycerol backbone

        # If all conditions are met
        return True, "Contains glycerol backbone esterified with two fatty acids and one phosphate group"

    return False, "Does not match phosphatidic acid structure"