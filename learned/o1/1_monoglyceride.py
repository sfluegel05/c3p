"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI:17854 1-monoglyceride
"""

from rdkit import Chem

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a monoglyceride in which the acyl substituent is located at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure there is exactly one ester group in the molecule
    ester_smarts = '[CX3](=O)[OX2H0][#6]'
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, requires exactly 1"

    # Get the ester oxygen atom
    ester_oxygen_idx = ester_matches[0][2]
    ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)

    # Traverse from the ester oxygen to find the glycerol backbone
    visited_atoms = set()
    atoms_to_visit = [(ester_oxygen, None)]  # (current_atom, previous_atom)

    glycerol_carbons = []
    hydroxyl_oxygens = []
    while atoms_to_visit:
        current_atom, previous_atom = atoms_to_visit.pop()
        if current_atom.GetIdx() in visited_atoms:
            continue
        visited_atoms.add(current_atom.GetIdx())

        # Check for carbon atoms connected to the ester oxygen
        if current_atom.GetAtomicNum() == 6:  # Carbon atom
            glycerol_carbons.append(current_atom)
            neighbors = current_atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetIdx() == previous_atom.GetIdx() if previous_atom else -1:
                    continue
                if neighbor.GetAtomicNum() == 6 or neighbor.GetAtomicNum() == 8:
                    atoms_to_visit.append((neighbor, current_atom))
        elif current_atom.GetAtomicNum() == 8:  # Oxygen atom
            # Check if it's a hydroxyl group
            if current_atom.GetTotalNumHs() == 1:
                hydroxyl_oxygens.append(current_atom)
                neighbors = current_atom.GetNeighbors()
                for neighbor in neighbors:
                    if neighbor.GetIdx() == previous_atom.GetIdx() if previous_atom else -1:
                        continue
                    if neighbor.GetAtomicNum() == 6:
                        atoms_to_visit.append((neighbor, current_atom))
            else:
                # Non-hydroxyl oxygen (e.g., part of ester), ignore
                continue

    # Check if glycerol backbone has exactly 3 carbons
    if len(glycerol_carbons) != 3:
        return False, "Glycerol backbone not identified correctly"

    # Check for correct substitutions on glycerol carbons
    # Sort carbons based on their distance from ester oxygen
    ester_carbon_idx = ester_matches[0][0]
    ester_carbon = mol.GetAtomWithIdx(ester_carbon_idx)

    # Map from atom index to distance from ester oxygen
    distances = Chem.GetDistanceMatrix(mol)
    distances_from_ester_oxygen = distances[ester_oxygen_idx]
    glycerol_carbons.sort(key=lambda atom: distances_from_ester_oxygen[atom.GetIdx()])

    # Position 1 carbon (closest to ester oxygen)
    pos1_carbon = glycerol_carbons[0]
    if pos1_carbon.GetDegree() != 2:
        return False, "Position 1 carbon has incorrect number of substituents"

    # Position 2 carbon
    pos2_carbon = glycerol_carbons[1]
    pos2_expected_substituents = set([pos1_carbon.GetIdx(), glycerol_carbons[2].GetIdx()])
    pos2_actual_substituents = set([atom.GetIdx() for atom in pos2_carbon.GetNeighbors()])
    if not pos2_expected_substituents.issubset(pos2_actual_substituents):
        return False, "Position 2 carbon not connected to positions 1 and 3 carbons"

    # Check that position 2 carbon has hydroxyl group
    pos2_has_oh = False
    for neighbor in pos2_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
            pos2_has_oh = True
    if not pos2_has_oh:
        return False, "Position 2 carbon does not have hydroxyl group"

    # Position 3 carbon
    pos3_carbon = glycerol_carbons[2]
    if pos3_carbon.GetDegree() != 2:
        return False, "Position 3 carbon has incorrect number of substituents"

    # Check that position 3 carbon has hydroxyl group
    pos3_has_oh = False
    for neighbor in pos3_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:
            pos3_has_oh = True
    if not pos3_has_oh:
        return False, "Position 3 carbon does not have hydroxyl group"

    # Check that there are no additional substituents on glycerol carbons
    for carbon in glycerol_carbons:
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() not in [6, 8]:  # Carbon or oxygen
                return False, "Additional substituents found on glycerol backbone"

    # Ensure there are no phosphate or amine groups attached to the molecule
    # Look for phosphate groups
    phosphate_smarts = '[P](=O)([O-])[O]'
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group found in molecule"

    # Look for amine groups
    amine_smarts = '[NX3;H2,H1;!$(NC=O)]'
    amine_pattern = Chem.MolFromSmarts(amine_smarts)
    if mol.HasSubstructMatch(amine_pattern):
        return False, "Amine group found in molecule"

    return True, "Contains glycerol backbone with acyl group at position 1"