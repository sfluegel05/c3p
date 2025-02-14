"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:64675 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is a lipid similar in structure to a glycerolipid but in which
    one or more of the carbon atoms on glycerol is bonded to an alkyl chain via an
    ether linkage, as opposed to the usual ester linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Glycerol backbone: O-C-C-C-O with hydroxyl groups
    glycerol_pattern = Chem.MolFromSmarts("""
    [O;!R][CH2][CH](O)[CH2][O;!R]
    """)

    # Ether linkage: C-O-C where O is connected to glycerol carbon
    ether_linkage_pattern = Chem.MolFromSmarts("""
    [C;!R][O;!R][CH2][CH](O)[CH2][O;!R]  # Ether linkage at primary carbon
    |
    [O;!R][CH2][CH](O)[CH2][O;!R][C;!R]  # Ether linkage at terminal carbon
    |
    [O;!R][CH2][CH](O)[CH2][O;!R]
    """)

    # Match glycerol backbone
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"

    # Initialize flags
    ether_linkage_found = False
    ester_linkage_found = False

    # Iterate over glycerol backbones found
    for match in glycerol_matches:
        # Get atoms in the glycerol backbone
        o1_idx, c1_idx, c2_idx, c3_idx, o2_idx = match
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)
        o1 = mol.GetAtomWithIdx(o1_idx)
        o2 = mol.GetAtomWithIdx(o2_idx)

        # Check for ether linkages on carbons c1, c2, c3
        # Function to check for ether linkage on a given carbon
        def has_ether_linkage(carbon_atom):
            for nbr in carbon_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in [o1_idx, o2_idx]:
                    # Check if oxygen is connected to another carbon (alkyl chain)
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetAtomicNum() == 6 and nbr2.GetIdx() != carbon_atom.GetIdx():
                            # Ensure it's an ether linkage (not ester)
                            bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nbr2.GetIdx())
                            if bond.GetBondType() == Chem.BondType.SINGLE:
                                # Check that nbr2 is part of an alkyl chain
                                return True
            return False

        # Check for ester linkage (C=O) on carbons c1, c2, c3
        def has_ester_linkage(carbon_atom):
            for nbr in carbon_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        return True
            return False

        # Check each carbon for ether or ester linkages
        c_atoms = [c1, c2, c3]
        for c in c_atoms:
            if has_ether_linkage(c):
                ether_linkage_found = True
            if has_ester_linkage(c):
                ester_linkage_found = True

    if not ether_linkage_found:
        return False, "No ether linkage found on glycerol backbone"

    # Since it's similar in structure to a glycerolipid, it can have ester linkages too
    return True, "Contains glycerol backbone with at least one alkyl chain attached via ether linkage"