"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: monoacylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string.
    A monoacylglycerol is a glyceride in which any one of the R groups (position not specified)
    is an acyl group while the remaining two R groups can be either H or alkyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacylglycerol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone pattern: three carbons in a chain, each connected to an oxygen
    glycerol_pattern = Chem.MolFromSmarts("[CH2](-[O])[CH](-[O])[CH2](-[O])")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone with oxygens found"

    # For each glycerol backbone found
    for match in matches:
        # Atoms in the match: (C1, O1, C2, O2, C3, O3)
        C1 = mol.GetAtomWithIdx(match[0])
        O1 = mol.GetAtomWithIdx(match[1])
        C2 = mol.GetAtomWithIdx(match[2])
        O2 = mol.GetAtomWithIdx(match[3])
        C3 = mol.GetAtomWithIdx(match[4])
        O3 = mol.GetAtomWithIdx(match[5])

        ester_count = 0
        other_oxygens = []

        # Function to check if an oxygen is part of an ester linkage
        def is_ester_oxygen(oxygen_atom, glycerol_carbon_idx):
            for bond in oxygen_atom.GetBonds():
                neighbor = bond.GetOtherAtom(oxygen_atom)
                if neighbor.GetIdx() == glycerol_carbon_idx:
                    continue  # Skip the glycerol carbon
                if neighbor.GetAtomicNum() == 6:  # Carbon
                    # Check if the carbon is a carbonyl carbon (C=O)
                    has_carbonyl = False
                    for nb in neighbor.GetNeighbors():
                        if nb.GetAtomicNum() == 8:
                            bond_type = neighbor.GetBondBetweenAtoms(neighbor.GetIdx(), nb.GetIdx()).GetBondType()
                            if bond_type == Chem.BondType.DOUBLE:
                                has_carbonyl = True
                                break
                    if has_carbonyl:
                        return True
            return False

        # Check O1
        if is_ester_oxygen(O1, C1.GetIdx()):
            ester_count += 1
        else:
            other_oxygens.append((O1, C1.GetIdx()))

        # Check O2
        if is_ester_oxygen(O2, C2.GetIdx()):
            ester_count += 1
        else:
            other_oxygens.append((O2, C2.GetIdx()))

        # Check O3
        if is_ester_oxygen(O3, C3.GetIdx()):
            ester_count += 1
        else:
            other_oxygens.append((O3, C3.GetIdx()))

        if ester_count != 1:
            continue  # Not a monoacylglycerol, try next match

        # Check that the other two oxygens are connected to H or alkyl groups
        valid_substituents = True
        for oxygen_atom, glycerol_carbon_idx in other_oxygens:
            neighbors = oxygen_atom.GetNeighbors()
            substituent_count = 0
            for neighbor in neighbors:
                if neighbor.GetIdx() == glycerol_carbon_idx:
                    continue  # Skip the glycerol carbon
                elif neighbor.GetAtomicNum() == 1:
                    # Hydrogen atom (OH group)
                    continue
                elif neighbor.GetAtomicNum() == 6:
                    # Carbon atom (alkyl group)
                    # Check that it's not a carbonyl carbon (not ester or acid)
                    is_carbonyl = False
                    for nb in neighbor.GetNeighbors():
                        if nb.GetAtomicNum() == 8:
                            bond_type = neighbor.GetBondBetweenAtoms(neighbor.GetIdx(), nb.GetIdx()).GetBondType()
                            if bond_type == Chem.BondType.DOUBLE:
                                is_carbonyl = True
                                break
                    if is_carbonyl:
                        valid_substituents = False
                        break
                else:
                    # Other atoms are not allowed
                    valid_substituents = False
                    break
            if not valid_substituents:
                break

        if not valid_substituents:
            continue  # Not a monoacylglycerol, try next match

        # Passed all checks
        return True, "Contains glycerol backbone with one acyl group and two H/alkyl groups"

    # If none of the matches meet the criteria
    return False, "Does not match monoacylglycerol pattern"