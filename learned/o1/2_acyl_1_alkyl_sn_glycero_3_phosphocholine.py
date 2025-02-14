"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: CHEBI:63866 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.
    This class has:
    - A glycerol backbone with sn (stereospecific numbering) configuration at position 2.
    - An ether-linked alkyl chain at position 1.
    - An ester-linked acyl chain at position 2.
    - A phosphocholine group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for the required substructures with atom mapping

    # Glycerol backbone with sn stereochemistry at position 2
    # Atoms are labeled with mapping numbers for positions
    glycerol_sn_pattern = Chem.MolFromSmarts("[C:1][C@H:2](O)[C:3]")
    if not mol.HasSubstructMatch(glycerol_sn_pattern):
        return False, "No glycerol backbone with sn stereochemistry found"

    # Ether-linked alkyl chain at position 1
    # Position 1 carbon connected to oxygen and alkyl chain
    ether_alkyl_pattern = Chem.MolFromSmarts("[C:1]-[O]-[C;!$(C=O)]")
    if not mol.HasSubstructMatch(ether_alkyl_pattern):
        return False, "No ether-linked alkyl chain at position 1 found"

    # Ester-linked acyl chain at position 2
    # Position 2 oxygen connected to carbonyl carbon
    ester_acyl_pattern = Chem.MolFromSmarts("[C@H:2]-[O]-[C:4](=O)-[C]")
    if not mol.HasSubstructMatch(ester_acyl_pattern):
        return False, "No ester-linked acyl chain at position 2 found"

    # Phosphocholine group at position 3
    phosphocholine_pattern = Chem.MolFromSmarts("[C:3]-O-P(=O)([O-])-OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group at position 3 found"

    # Ensure that the patterns are matched at the correct positions
    # Get the matches for glycerol backbone
    backbone_matches = mol.GetSubstructMatches(glycerol_sn_pattern)
    for match in backbone_matches:
        c1_idx = match[0]
        c2_idx = match[1]
        c3_idx = match[2]

        # Check that ether-linked alkyl chain is connected to position 1 carbon
        c1_atom = mol.GetAtomWithIdx(c1_idx)
        ether_alkyl_match = False
        for bond in c1_atom.GetBonds():
            nbr = bond.GetOtherAtom(c1_atom)
            if nbr.GetAtomicNum() == 8:  # Oxygen
                # Check if oxygen is connected to an alkyl group
                for obond in nbr.GetBonds():
                    onbr = obond.GetOtherAtom(nbr)
                    if onbr.GetAtomicNum() == 6 and onbr.GetIdx() != c1_idx:
                        ether_alkyl_match = True
                        break
            if ether_alkyl_match:
                break
        if not ether_alkyl_match:
            continue  # Try next backbone match

        # Check that ester-linked acyl chain is connected to position 2 oxygen
        c2_atom = mol.GetAtomWithIdx(c2_idx)
        ester_acyl_match = False
        for bond in c2_atom.GetBonds():
            nbr = bond.GetOtherAtom(c2_atom)
            if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != c3_idx:  # Oxygen not connected to position 3 carbon
                # Oxygen connected to carbonyl carbon
                for obond in nbr.GetBonds():
                    onbr = obond.GetOtherAtom(nbr)
                    if onbr.GetAtomicNum() == 6 and obond.GetBondTypeAsDouble() == 1.0:
                        # Check if carbon is carbonyl carbon
                        has_double_bonded_oxygen = False
                        for cbond in onbr.GetBonds():
                            cnbr = cbond.GetOtherAtom(onbr)
                            if cnbr.GetAtomicNum() == 8 and cbond.GetBondTypeAsDouble() == 2.0:
                                has_double_bonded_oxygen = True
                                break
                        if has_double_bonded_oxygen:
                            ester_acyl_match = True
                            break
            if ester_acyl_match:
                break
        if not ester_acyl_match:
            continue  # Try next backbone match

        # Check that phosphocholine group is connected to position 3 carbon
        c3_atom = mol.GetAtomWithIdx(c3_idx)
        phosphocholine_match = False
        for bond in c3_atom.GetBonds():
            nbr = bond.GetOtherAtom(c3_atom)
            if nbr.GetAtomicNum() == 8:  # Oxygen
                # Oxygen connected to phosphate
                for obond in nbr.GetBonds():
                    onbr = obond.GetOtherAtom(nbr)
                    if onbr.GetAtomicNum() == 15:  # Phosphorus
                        phosphocholine_match = True
                        break
            if phosphocholine_match:
                break
        if not phosphocholine_match:
            continue  # Try next backbone match

        # All checks passed
        return True, "Molecule matches all required substructures for 2-acyl-1-alkyl-sn-glycero-3-phosphocholine"

    return False, "Molecule does not match all required substructures"