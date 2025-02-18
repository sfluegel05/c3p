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

    # Check for presence of an amino group (NH2 or NH3+)
    amino_groups = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Check for NH2 or NH3+
            if (atom.GetFormalCharge() == 0 and atom.GetTotalNumHs() >= 2) or \
               (atom.GetFormalCharge() == 1 and atom.GetTotalNumHs() >= 3):
                amino_groups.append(atom)
    if not amino_groups:
        return False, "No amino group (NH2/NH3+) found"

    # Check if any amino group has a hydroxyl group on adjacent carbon
    has_adjacent_oh = False
    for amino_atom in amino_groups:
        for neighbor in amino_atom.GetNeighbors():
            # Check neighbor's bonds for hydroxyl (O with at least one H)
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

    # Check for long carbon chain (at least 12 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 12:
        return False, f"Carbon count ({carbon_count}) too low for sphingoid"

    # Check for at least one additional hydroxyl group (1,3-diol or similar)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1)
    if hydroxyl_count < 2:
        return False, "Insufficient hydroxyl groups for sphingoid"

    return True, "Long-chain amino alcohol with adjacent hydroxyl groups"