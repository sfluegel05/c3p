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

    # Define SMARTS patterns for the glycerol backbone
    # Glycerol backbone with hydroxyls at positions 2 and 3
    glycerol_smarts = '[C@H](O)[C@@H](O)CO'  # Represents glycerol backbone
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)

    if not glycerol_matches:
        return False, "No glycerol backbone with hydroxyls at positions 2 and 3 found"

    # Define SMARTS pattern for ester group attached to primary carbon (position 1)
    ester_smarts = 'O[C@@H](CO)COC(=O)[C;!H0]'  # Ester linked at position 1
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    if not ester_matches:
        return False, "No ester group attached at position 1 found"

    # Check for only one ester group
    ester_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.ESTER:
            ester_count +=1
    if ester_count > 1:
        return False, f"More than one ester group found ({ester_count} esters)"

    # Ensure that the ester group is at position 1
    # Find the atom indices for glycerol carbons and ester oxygen
    for match in glycerol_matches:
        # Indices for glycerol backbone carbons
        c2_idx = match[0]
        c3_idx = match[1]
        c1_idx = match[2]

        c1_atom = mol.GetAtomWithIdx(c1_idx)
        c2_atom = mol.GetAtomWithIdx(c2_idx)
        c3_atom = mol.GetAtomWithIdx(c3_idx)

        # Check if c1 is connected via an ester linkage
        is_esterified = False
        for neighbor in c1_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for nbr in neighbor.GetNeighbors():
                    if nbr.GetIdx() != c1_idx and nbr.GetAtomicNum() == 6:
                        # Check if it's a carbonyl carbon
                        for bond in neighbor.GetBonds():
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetOtherAtomIdx(neighbor.GetIdx()) == nbr.GetIdx():
                                is_esterified = True
                                break
        if is_esterified:
            # Verify that c2 and c3 have hydroxyl groups
            c2_oh = False
            c3_oh = False
            for neighbor in c2_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    c2_oh = True
                    break
            for neighbor in c3_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                    c3_oh = True
                    break
            if c2_oh and c3_oh:
                return True, "Contains glycerol backbone with acyl group at position 1"
            else:
                return False, "Positions 2 and/or 3 do not have hydroxyl groups"

    return False, "No esterified glycerol backbone at position 1 found"