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

    # Define ester pattern: O=C-O-C 
    ester_pattern = Chem.MolFromSmarts('C(=O)O[*]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Collect the indices of ester oxygens and the carbons they are attached to (glycerol carbons)
    esterified_carbons = set()
    for match in ester_matches:
        carbonyl_c = match[0]   # Carbonyl carbon
        ester_o = match[1]      # Ester oxygen
        attached_atom = match[2]  # Atom attached to ester oxygen (should be glycerol carbon)

        # Check if the attached atom is a carbon
        atom = mol.GetAtomWithIdx(attached_atom)
        if atom.GetAtomicNum() == 6:
            esterified_carbons.add(attached_atom)

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
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                # Check if this oxygen is part of an ester group
                bond = mol.GetBondBetweenAtoms(c_idx, neighbor.GetIdx())
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    for nbr in neighbor.GetNeighbors():
                        if nbr.GetIdx() != c_idx and nbr.GetAtomicNum() == 6:
                            # Found C-O-C where O is connected to a carbonyl carbon
                            for nbr_of_nbr in nbr.GetNeighbors():
                                if nbr_of_nbr.GetIdx() != neighbor.GetIdx() and nbr_of_nbr.GetAtomicNum() == 8:
                                    # Found carbonyl oxygen
                                    ester_o_counts +=1
                                    break
        if ester_o_counts != 1:
            return False, f"Glycerol carbon at index {c_idx} is not connected to exactly one ester group"

    # If all checks pass, the molecule is a triglyceride
    return True, "Molecule has a glycerol backbone with three fatty acid chains attached via ester bonds"