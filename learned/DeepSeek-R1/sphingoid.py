"""
Classifies: CHEBI:35785 sphingoid
"""
"""
Classifies: CHEBI: ??? sphingoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sphingoid(smiles: str):
    """
    Determines if a molecule is a sphingoid based on its SMILES string.
    Sphingoids are characterized by a long-chain amino alcohol backbone with hydroxyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of an amino group (NH2, NH3+, or amide)
    amino_groups = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Check for NH2, NH3+, or amide (N connected to carbonyl)
            if (atom.GetFormalCharge() in (0, 1) and atom.GetTotalNumHs() >= 1) or \
               any(bond.GetBondType() == Chem.BondType.SINGLE and bond.GetOtherAtom(atom).GetAtomicNum() == 6 and
                   any(nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 0 for nbr in bond.GetOtherAtom(atom).GetNeighbors() if nbr.GetAtomicNum() == 8)
                   for bond in atom.GetBonds()):
                amino_groups.append(atom)
    if not amino_groups:
        return False, "No amino group found"

    # Check for at least two hydroxyl groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 2:
        return False, "Insufficient hydroxyl groups (need at least two)"

    # Check for hydroxyl adjacent to any amino group
    has_adjacent_oh = False
    for amino_atom in amino_groups:
        for neighbor in amino_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                # Check if this carbon has a hydroxyl group
                for bond in neighbor.GetBonds():
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetAtomicNum() == 8 and other_atom.GetTotalNumHs() >= 1:
                        has_adjacent_oh = True
                        break
                if has_adjacent_oh:
                    break
        if has_adjacent_oh:
            break
    if not has_adjacent_oh:
        return False, "No hydroxyl group adjacent to amino group"

    # Check for long carbon chain (at least 12 carbons in the backbone)
    # Find the longest carbon chain
    chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            # Follow the chain in both directions
            current = atom
            visited = set()
            length = 1
            # Move in one direction
            while True:
                next_atoms = [nbr for nbr in current.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr not in visited]
                if len(next_atoms) != 1:
                    break
                visited.add(current)
                current = next_atoms[0]
                length += 1
            # Move in the other direction
            current = atom
            visited = set()
            while True:
                next_atoms = [nbr for nbr in current.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr not in visited]
                if len(next_atoms) != 1:
                    break
                visited.add(current)
                current = next_atoms[0]
                length += 1
            chain_length = max(chain_length, length)
    if chain_length < 12:
        return False, f"Longest carbon chain ({chain_length}) too short for sphingoid"

    return True, "Long-chain amino alcohol with adjacent hydroxyl groups"