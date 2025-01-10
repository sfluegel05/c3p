"""
Classifies: CHEBI:35759 1-monoglyceride
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone pattern with free hydroxyls at positions 2 and 3
    glycerol_pattern = Chem.MolFromSmarts("[CH2][CH](O)[CH2]O")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if len(glycerol_matches) == 0:
        return False, "Glycerol backbone not found"

    # Define ester group pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check if ester group is attached to position 1 of glycerol
    for glycerol_match in glycerol_matches:
        c1_idx = glycerol_match[0]
        c2_idx = glycerol_match[1]
        c3_idx = glycerol_match[2]
        o3_idx = glycerol_match[3]  # Oxygen attached to C3

        c1_atom = mol.GetAtomWithIdx(c1_idx)
        c2_atom = mol.GetAtomWithIdx(c2_idx)
        c3_atom = mol.GetAtomWithIdx(c3_idx)
        o3_atom = mol.GetAtomWithIdx(o3_idx)

        # Check that C2 and C3 have hydroxyl groups
        c2_has_oh = False
        c3_has_oh = False

        for neighbor in c2_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                bond = mol.GetBondBetweenAtoms(c2_atom.GetIdx(), neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE and neighbor.GetDegree() == 1:
                    c2_has_oh = True

        for neighbor in c3_atom.GetNeighbors():
            if neighbor.GetIdx() == o3_atom.GetIdx():
                # Confirm that O3 is a hydroxyl group
                if o3_atom.GetDegree() == 1:
                    c3_has_oh = True

        if not c2_has_oh or not c3_has_oh:
            continue  # This glycerol backbone doesn't have OH groups at positions 2 and 3

        # Check if C1 is connected to ester oxygen
        c1_has_ester = False
        for neighbor in c1_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                o_atom = neighbor
                # Check if this oxygen is part of the ester group
                for ester_match in ester_matches:
                    if o_atom.GetIdx() in ester_match:
                        c1_has_ester = True
                        break
            if c1_has_ester:
                break

        if c1_has_ester:
            return True, "Molecule is a 1-monoglyceride with acyl group at position 1 and free hydroxyls at positions 2 and 3"

    return False, "Molecule does not match 1-monoglyceride structure"