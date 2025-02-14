"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate
    based on its SMILES string.
    The molecule should have a glycerol backbone with a phosphate group at position 3,
    one acyl chain attached via an ester linkage at position 1 or 2,
    and no additional substitutions on the phosphate group.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule matches the criteria, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the glycerol phosphate backbone pattern (ignore stereochemistry)
    # Glycerol backbone with phosphate group at position 3
    glycerol_phosphate_smarts_list = [
        'OCC(O)COP(=O)(O)O',         # Protonated phosphate
        'OCC(O)COP(=O)([O-])O',      # Mono deprotonated phosphate
        'OCC(O)COP(=O)([O-])[O-]',   # Fully deprotonated phosphate
        'OCC(O)COP(=O)(O)[O-]',      # Mixed protonation
    ]

    # Try matching the glycerol phosphate backbone with different SMARTS patterns
    glycerol_phosphate_pattern = None
    for smarts in glycerol_phosphate_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        if matches:
            glycerol_phosphate_pattern = pattern
            break

    if glycerol_phosphate_pattern is None or not matches:
        return False, "No glycerol phosphate backbone found"

    # For each match, check for acyl chain and phosphate substitutions
    for match in matches:
        # Indices for glycerol carbons
        c1_idx = match[1]  # Carbon attached to primary hydroxyl or ester
        c2_idx = match[2]  # Central carbon with hydroxyl or ester
        c3_idx = match[4]  # Carbon attached to phosphate group

        # Get the phosphate group atom indices
        p_atom = mol.GetAtomWithIdx(match[5])
        if p_atom.GetAtomicNum() != 15:
            continue  # Not a phosphate group

        # Ensure phosphate has no substitutions other than oxygens
        phosphate_neighbors = [nbr.GetAtomicNum() for nbr in p_atom.GetNeighbors()]
        if any(num not in [8] for num in phosphate_neighbors):
            continue  # Phosphate is substituted

        # Check acyl chains at positions 1 and 2
        acyl_positions = []
        hydroxyl_positions = []
        for idx in [c1_idx, c2_idx]:
            carbon = mol.GetAtomWithIdx(idx)
            oxygen_atoms = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() == 8]
            ester_found = False
            hydroxyl_found = False
            for oxygen in oxygen_atoms:
                bonds = oxygen.GetBonds()
                double_bonded_c = None
                for bond in bonds:
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other_atom = bond.GetOtherAtom(oxygen)
                        if other_atom.GetAtomicNum() == 6:
                            double_bonded_c = other_atom
                            break
                if double_bonded_c:
                    # Found ester linkage
                    ester_found = True
                    break
                elif oxygen.GetTotalNumHs() > 0:
                    # Found hydroxyl group
                    hydroxyl_found = True
                    break
            if ester_found:
                acyl_positions.append(idx)
            elif hydroxyl_found:
                hydroxyl_positions.append(idx)
            else:
                # No ester or hydroxyl found at this position
                continue

        # Valid if we have exactly one acyl chain and one hydroxyl group at positions 1 and 2
        if len(acyl_positions) == 1 and len(hydroxyl_positions) == 1:
            # Check that acyl chain is sufficiently long (e.g., more than 3 carbons)
            acyl_carbon = mol.GetAtomWithIdx(acyl_positions[0])
            ester_oxygen = [nbr for nbr in acyl_carbon.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in match][0]
            next_atom = [nbr for nbr in ester_oxygen.GetNeighbors() if nbr.GetIdx() != acyl_carbon.GetIdx()][0]
            # Traverse the acyl chain
            acyl_length = 0
            visited = set()
            stack = [next_atom]
            while stack:
                current_atom = stack.pop()
                if current_atom.GetIdx() in visited:
                    continue
                visited.add(current_atom.GetIdx())
                if current_atom.GetAtomicNum() == 6:
                    acyl_length += 1
                    neighbors = [nbr for nbr in current_atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited]
                    stack.extend(neighbors)
            if acyl_length >= 4:
                return True, "Molecule matches monoacyl-sn-glycerol 3-phosphate structure"
            else:
                continue  # Acyl chain too short

    return False, "Molecule does not match the required structure"