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
    
    # Define the glycerol phosphate backbone pattern
    # Glycerol backbone with phosphate group at position 3
    glycerol_phosphate_pattern = Chem.MolFromSmarts('[C@@H](O)[C@@H](O)COP(O)(O)=O')
    if glycerol_phosphate_pattern is None:
        return False, "Invalid SMARTS pattern for glycerol phosphate"

    matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
    if not matches:
        # Try matching without stereochemistry
        glycerol_phosphate_pattern = Chem.MolFromSmarts('C(O)C(O)COP(O)(O)=O')
        matches = mol.GetSubstructMatches(glycerol_phosphate_pattern)
        if not matches:
            return False, "No glycerol phosphate backbone found"

    # For each match, check for acyl chain and phosphate substitutions
    for match in matches:
        c1_idx = match[0]
        c2_idx = match[1]
        c3_idx = match[2]
        o_p_idx = match[3]  # Oxygen connecting to phosphate

        # Ensure phosphate group is not substituted
        phosphate_atom = mol.GetAtomWithIdx(o_p_idx).GetNeighbors()[0]  # Get the phosphorus atom
        if phosphate_atom.GetAtomicNum() != 15:
            return False, "Phosphate group not found"
        phosphate_neighbors = phosphate_atom.GetNeighbors()
        non_oxygen_neighbors = [a for a in phosphate_neighbors if a.GetAtomicNum() != 8 and a.GetIdx() != o_p_idx]
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
                if neighbor.GetAtomicNum() == 8:
                    # Check if oxygen is part of an ester linkage
                    oxygen = neighbor
                    for obond in oxygen.GetBonds():
                        o_neighbor = obond.GetOtherAtom(oxygen)
                        if obond.GetBondType() == Chem.rdchem.BondType.DOUBLE and o_neighbor.GetAtomicNum() == 6:
                            # Ester linkage found
                            ester_found = True
                            break
                        elif obond.GetBondType() == Chem.rdchem.BondType.SINGLE and o_neighbor.GetAtomicNum() == 1:
                            # Hydroxyl group
                            hydroxyl_found = True
                    if ester_found or hydroxyl_found:
                        break
            if ester_found:
                acyl_positions += 1
            elif hydroxyl_found:
                hydroxyl_positions +=1
            else:
                # Neither ester nor hydroxyl found
                continue

        # Valid if we have exactly one acyl chain and one hydroxyl group at positions 1 and 2
        if acyl_positions == 1 and hydroxyl_positions ==1:
            return True, "Molecule matches monoacyl-sn-glycerol 3-phosphate structure"
    
    return False, "Molecule does not match the required structure"