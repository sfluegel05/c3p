"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate or diacyl-sn-glycerol 3-phosphate
    based on its SMILES string.
    The molecule should have a glycerol backbone with a phosphate group at position 3,
    one or two acyl chains attached via ester linkages at positions 1 and/or 2,
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
    
    # Define the glycerol phosphate backbone pattern
    # Glycerol backbone with phosphate group at position 3
    glycerol_phosphate_pattern = Chem.MolFromSmarts("""
    [C@@H](O)[C@@H](O)[C@@H](O[P](=O)(O)O)
    """)
    if glycerol_phosphate_pattern is None:
        return False, "Invalid SMARTS pattern for glycerol phosphate"

    matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    if not matches:
        # Try matching without stereochemistry
        glycerol_phosphate_pattern = Chem.MolFromSmarts("""
        [C](O)[C](O)[C](O[P](=O)(O)O)
        """)
        matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
        if not matches:
            return False, "No glycerol phosphate backbone found"

    # For each match, check for acyl chains and phosphate substitutions
    for match in matches:
        c1_idx = match[0]
        c2_idx = match[1]
        c3_idx = match[2]
        phosphate_o_idx = match[3]
        phosphate_p_idx = match[4]

        # Ensure phosphate group is not substituted
        phosphate_atom = mol.GetAtomWithIdx(phosphate_p_idx)
        phosphate_neighbors = phosphate_atom.GetNeighbors()
        non_oxygen_neighbors = [a for a in phosphate_neighbors if a.GetAtomicNum() != 8]
        if len(non_oxygen_neighbors) > 0:
            # Phosphate is substituted
            continue

        # Check acyl chains at positions 1 and 2
        acyl_positions = 0
        hydroxyl_positions = 0
        for idx in [c1_idx, c2_idx]:
            carbon = mol.GetAtomWithIdx(idx)
            ester_found = False
            hydroxyl_found = False
            for bond in carbon.GetBonds():
                neighbor = bond.GetOtherAtom(carbon)
                if neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Check if oxygen is part of an ester or hydroxyl
                    oxygen = neighbor
                    connected_atoms = [a.GetAtomicNum() for a in oxygen.GetNeighbors()]
                    if connected_atoms.count(6) > 1:
                        # Part of an ether linkage, not expected here
                        continue
                    elif any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and
                             a.GetAtomicNum() == 6 for a in oxygen.GetNeighbors() for bond in a.GetBonds()):
                        # Part of an ester linkage
                        ester_found = True
                        break
                    elif len(oxygen.GetNeighbors()) == 1:
                        # Hydroxyl group
                        hydroxyl_found = True
            if ester_found:
                acyl_positions += 1
            elif hydroxyl_found:
                hydroxyl_positions +=1

        # Valid if we have one or two acyl chains and remaining positions are hydroxyls
        if acyl_positions >= 1 and acyl_positions + hydroxyl_positions == 2:
            return True, "Molecule matches monoacyl-sn-glycerol 3-phosphate structure"
    
    return False, "Molecule does not match the required structure"