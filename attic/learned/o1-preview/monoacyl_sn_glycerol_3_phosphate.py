"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    A monoacyl-sn-glycerol 3-phosphate is a glycerol backbone with a phosphate group at position 3,
    and a single acyl group attached via ester linkage at either position 1 or position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule is sanitized
    try:
        Chem.SanitizeMol(mol)
    except Chem.rdchem.KekulizeException:
        return False, "Molecule could not be sanitized"

    # Define phosphate group pattern connected to carbon via oxygen
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)(O)O")
    if phosphate_pattern is None:
        return False, "Invalid phosphate SMARTS pattern"

    # Find phosphate groups connected to glycerol backbone
    matches = mol.GetSubstructMatches(phosphate_pattern)
    if not matches:
        return False, "No phosphate group found"

    # Initialize variables
    num_acyl_groups = 0
    glycerol_carbons = set()

    for match in matches:
        # Get phosphate atom index
        phosphate_idx = match[0]
        phosphate_atom = mol.GetAtomWithIdx(phosphate_idx)

        # Find the oxygen connecting phosphate to glycerol carbon (position 3)
        phosphate_oxygen = None
        for neighbor in phosphate_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:
                # Check if this oxygen is connected to a carbon (glycerol C3)
                for nbr in neighbor.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != phosphate_atom.GetIdx():
                        phosphate_oxygen = neighbor
                        glycerol_c3 = nbr
                        break
                if phosphate_oxygen:
                    break
        
        if not phosphate_oxygen:
            continue  # No glycerol backbone connected to phosphate

        # Collect glycerol carbons (C1, C2, C3)
        glycerol_c3_idx = glycerol_c3.GetIdx()
        glycerol_carbons.add(glycerol_c3_idx)

        # Get C2 (connected to C3)
        c2_atoms = [nbr for nbr in glycerol_c3.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != phosphate_oxygen.GetIdx()]
        if not c2_atoms:
            continue  # No C2 found
        glycerol_c2 = c2_atoms[0]
        glycerol_carbons.add(glycerol_c2.GetIdx())

        # Get C1 (connected to C2)
        c1_atoms = [nbr for nbr in glycerol_c2.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != glycerol_c3_idx]
        if not c1_atoms:
            continue  # No C1 found
        glycerol_c1 = c1_atoms[0]
        glycerol_carbons.add(glycerol_c1.GetIdx())

        # Check substituents at C1 and C2
        acyl_positions = []
        hydroxyl_positions = []

        for carbon_atom in [glycerol_c1, glycerol_c2]:
            has_acyl = False
            has_oh = False
            for nbr in carbon_atom.GetNeighbors():
                # Skip bonds to other glycerol carbons
                if nbr.GetIdx() in glycerol_carbons:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Check if oxygen is part of ester linkage (C=O)-O-C
                    is_ester = False
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetAtomicNum() == 6 and nbr2.GetIdx() != carbon_atom.GetIdx():
                            for nbr3 in nbr2.GetNeighbors():
                                if nbr3.GetAtomicNum() == 8 and nbr3.GetIdx() != nbr.GetIdx() and mol.GetBondBetweenAtoms(nbr2.GetIdx(), nbr3.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    is_ester = True
                                    break
                    if is_ester:
                        has_acyl = True
                    else:
                        has_oh = True
                elif nbr.GetAtomicNum() == 6:
                    # Check for ethers or other linkages (ignore for now)
                    continue
            if has_acyl:
                acyl_positions.append(carbon_atom.GetIdx())
                num_acyl_groups += 1
            elif has_oh:
                hydroxyl_positions.append(carbon_atom.GetIdx())

        # Monoacyl-sn-glycerol 3-phosphate should have exactly one acyl group at position 1 or 2
        if num_acyl_groups != 1:
            return False, f"Expected exactly 1 acyl group at position 1 or 2, found {num_acyl_groups}"

        if len(hydroxyl_positions) != 1:
            return False, f"Expected exactly 1 hydroxyl group at position 1 or 2, found {len(hydroxyl_positions)}"

        return True, "Molecule is a monoacyl-sn-glycerol 3-phosphate"

    return False, "No glycerol 3-phosphate backbone found"