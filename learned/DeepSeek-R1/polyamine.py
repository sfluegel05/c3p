"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:53371 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine contains two or more amino groups (-NH2, -NH-, etc.) in an organic compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    amine_count = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Only process nitrogen atoms

        # Check for nitro groups (N with two double-bonded oxygens)
        nitro_oxygens = 0
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 8:
                    nitro_oxygens += 1
        if nitro_oxygens >= 2:
            continue

        # Check for nitriles (triple bond to carbon)
        is_nitrile = any(bond.GetBondType() == Chem.BondType.TRIPLE 
                        and bond.GetOtherAtom(atom).GetAtomicNum() == 6 
                        for bond in atom.GetBonds())
        if is_nitrile:
            continue

        # Check for amides (bonded to carbonyl group)
        is_amide = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                for bond in neighbor.GetBonds():
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(neighbor)
                        if other.GetAtomicNum() == 8:
                            is_amide = True
                            break
                if is_amide:
                    break
        if is_amide:
            continue

        # Check for sulfonamides (bonded to sulfonyl group)
        is_sulfonamide = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 16:  # Sulfur
                sulfur_bonds = neighbor.GetBonds()
                oxygen_double = sum(1 for b in sulfur_bonds 
                                   if b.GetBondType() == Chem.BondType.DOUBLE 
                                   and b.GetOtherAtom(neighbor).GetAtomicNum() == 8)
                if oxygen_double >= 2:
                    is_sulfonamide = True
                    break
        if is_sulfonamide:
            continue

        # Check if bonded to at least one carbon
        has_carbon = any(n.GetAtomicNum() == 6 for n in atom.GetNeighbors())
        if has_carbon:
            amine_count += 1

    if amine_count >= 2:
        return True, f"Contains {amine_count} amine groups"
    return False, f"Only {amine_count} amine groups found (requires ≥2)"