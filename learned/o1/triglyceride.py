"""
Classifies: CHEBI:17855 triglyceride
"""
"""
Classifies: CHEBI:17855 triglyceride
"""
from rdkit import Chem

def is_triglyceride(smiles: str):
    """
    Determines if a molecule is a triglyceride based on its SMILES string.
    A triglyceride is a glycerol backbone with three fatty acid chains attached via ester bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triglyceride, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define ester pattern: [#6][C](=O)[O][#6]
    ester_pattern = Chem.MolFromSmarts('[#6][CX3](=O)[OX2H0][#6]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Collect the indices of ester oxygens and the carbons they are attached to (glycerol carbons)
    esterified_carbons = set()
    for match in ester_matches:
        # match indices: [C(attached to carbonyl C), carbonyl C, ester O, glycerol C]
        glycerol_c = match[3]   # Carbon attached to ester oxygen (possibly glycerol carbon)
        esterified_carbons.add(glycerol_c)

    # Check if there are exactly 3 esterified carbons
    if len(esterified_carbons) != 3:
        return False, f"Found {len(esterified_carbons)} ester groups attached to carbons, need exactly 3"

    # List of esterified carbon indices
    glycerol_carbons = list(esterified_carbons)

    # Check connectivity between glycerol carbons
    bonds = []
    for i in range(3):
        for j in range(i+1, 3):
            bond = mol.GetBondBetweenAtoms(glycerol_carbons[i], glycerol_carbons[j])
            if bond:
                bonds.append((glycerol_carbons[i], glycerol_carbons[j]))

    # There should be exactly 2 bonds connecting the glycerol carbons (forming a chain)
    if len(bonds) != 2:
        return False, "Glycerol carbons are not connected properly"

    # Create a subgraph of glycerol carbons to analyze connectivity
    from collections import defaultdict
    connectivity = defaultdict(list)
    for a1, a2 in bonds:
        connectivity[a1].append(a2)
        connectivity[a2].append(a1)

    # The central carbon should be connected to two others
    central_carbons = [atom_idx for atom_idx, neighbors in connectivity.items() if len(neighbors) == 2]
    if len(central_carbons) != 1:
        return False, "Glycerol backbone does not have a central carbon connected to two others"

    # The other two carbons should each be connected to the central carbon only
    terminal_carbons = [atom_idx for atom_idx, neighbors in connectivity.items() if len(neighbors) == 1]
    if len(terminal_carbons) != 2:
        return False, "Glycerol backbone does not have two terminal carbons connected to the central carbon"

    # Check that each glycerol carbon is connected to exactly one ester oxygen
    for c_idx in glycerol_carbons:
        atom = mol.GetAtomWithIdx(c_idx)
        ester_o_counts = 0
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                # Check if this oxygen is part of an ester group
                o_atom = neighbor
                is_ester_o = False
                for o_bond in o_atom.GetBonds():
                    o_neighbor = o_bond.GetOtherAtom(o_atom)
                    if o_neighbor.GetIdx() != c_idx:
                        if o_neighbor.GetAtomicNum() == 6:
                            # Check if this carbon (o_neighbor) is a carbonyl carbon with a double-bonded oxygen
                            has_carbonyl_o = False
                            for obond in o_neighbor.GetBonds():
                                if obond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                    other_atom = obond.GetOtherAtom(o_neighbor)
                                    if other_atom.GetAtomicNum() == 8:
                                        has_carbonyl_o = True
                                        break
                            if has_carbonyl_o:
                                is_ester_o = True
                                break
                if is_ester_o:
                    ester_o_counts +=1
        if ester_o_counts !=1:
            return False, f"Glycerol carbon at index {c_idx} is not connected to exactly one ester group"

    # If all checks pass, the molecule is a triglyceride
    return True, "Molecule has a glycerol backbone with three fatty acid chains attached via ester bonds"