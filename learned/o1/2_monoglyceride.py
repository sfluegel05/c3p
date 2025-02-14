"""
Classifies: CHEBI:17389 2-monoglyceride
"""
"""
Classifies: 2-monoglyceride
"""

from rdkit import Chem

def is_2_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 2-monoglyceride based on its SMILES string.
    A 2-monoglyceride is a monoglyceride where the acyl group is attached at position 2 of glycerol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-monoglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol backbone pattern (three connected carbons)
    glycerol_pattern = Chem.MolFromSmarts("[C;!R]-[C;!R]-[C;!R]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)

    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Iterate over each glycerol backbone match
    for match in glycerol_matches:
        c1_idx, c2_idx, c3_idx = match

        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)

        # Check if c1 and c3 have hydroxyl groups (OH)
        def has_hydroxyl(carbon_atom, exclude_idx):
            for neighbor in carbon_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    # Oxygen with only one bond (hydroxyl)
                    if neighbor.GetIdx() != exclude_idx:
                        return True
            return False

        c1_has_oh = has_hydroxyl(c1, c2_idx)
        c3_has_oh = has_hydroxyl(c3, c2_idx)

        if not (c1_has_oh and c3_has_oh):
            continue  # Positions 1 and 3 must have hydroxyl groups

        # Check if c2 is connected to an ester group
        ester_found = False
        for neighbor in c2.GetNeighbors():
            if neighbor.GetIdx() == c1_idx or neighbor.GetIdx() == c3_idx:
                continue  # Skip the glycerol backbone bonds

            if neighbor.GetAtomicNum() == 8:
                # Check if this oxygen is part of an ester linkage
                for neighbor2 in neighbor.GetNeighbors():
                    if neighbor2.GetIdx() == c2_idx:
                        continue  # Skip the bond back to c2
                    if neighbor2.GetAtomicNum() == 6:
                        # Check for carbonyl (C=O)
                        for bond in neighbor2.GetBonds():
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetOtherAtomIdx(neighbor2.GetIdx()) != neighbor.GetIdx():
                                ester_found = True
                                break
                        if ester_found:
                            break
                if ester_found:
                    break
        if not ester_found:
            continue  # Position 2 must have an ester linkage

        # Check that there is only one ester group in the molecule
        ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
        ester_matches = mol.GetSubstructMatches(ester_pattern)
        if len(ester_matches) != 1:
            return False, f"Found {len(ester_matches)} ester groups, expected 1"

        return True, "Molecule is a 2-monoglyceride with acyl group at position 2"

    return False, "Molecule lacks required ester linkage at position 2 or hydroxyls at positions 1 and 3"