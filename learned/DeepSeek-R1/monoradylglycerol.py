"""
Classifies: CHEBI:76575 monoradylglycerol
"""
"""
Classifies: CHEBI monoradylglycerol (Any lipid that is glycerol bearing a single acyl, alkyl or alk-1-enyl substituent)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoradylglycerol(smiles: str):
    """
    Determines if a molecule is a monoradylglycerol based on its SMILES string.
    A monoradylglycerol has a glycerol backbone with a single acyl (ester), alkyl (ether), or alk-1-enyl (vinyl ether) substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoradylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # SMARTS pattern for glycerol backbone with two hydroxyls and one O-linked substituent (any position)
    glycerol_pattern = Chem.MolFromSmarts(
        "[C](-[OH])-[C](-[OH])-[C](-O[!H0]) | "  # substituent on third carbon
        "[C](-[OH])-[C](-O[!H0])-[C](-[OH]) | "  # substituent on second carbon
        "[C](-O[!H0])-[C](-[OH])-[C](-[OH])"     # substituent on first carbon
    )
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with two hydroxyls and one O-linked group"

    # Check for exactly one O-linked substituent (excluding hydroxyls)
    o_substituents = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            # Check if oxygen is part of a substituent (not hydroxyl)
            if not any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors()):
                # Check if connected to a carbon in the glycerol backbone
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.HasSubstructMatch(glycerol_pattern):
                        o_substituents += 1
                        break
    if o_substituents != 1:
        return False, f"Found {o_substituents} substituents, need exactly 1"

    # Check substituent type (acyl, alkyl, alk-1-enyl)
    matches = mol.GetSubstructMatches(glycerol_pattern)
    for match in matches:
        # Find the substituent oxygen
        substituent_oxygen = None
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and not any(n.GetAtomicNum() == 1 for n in neighbor.GetNeighbors()):
                    substituent_oxygen = neighbor
                    break
            if substituent_oxygen:
                break

        if not substituent_oxygen:
            continue

        # Check for ester (acyl)
        ester = False
        for bond in substituent_oxygen.GetBonds():
            if bond.GetBondType() == Chem.BondType.SINGLE:
                next_atom = bond.GetOtherAtom(substituent_oxygen)
                if next_atom.GetAtomicNum() == 6:
                    for next_bond in next_atom.GetBonds():
                        if next_bond.GetBondType() == Chem.BondType.DOUBLE:
                            other_atom = next_bond.GetOtherAtom(next_atom)
                            if other_atom.GetAtomicNum() == 8:
                                ester = True
                                break

        # Check for vinyl ether (alk-1-enyl)
        vinyl_ether = False
        if not ester:
            for bond in substituent_oxygen.GetBonds():
                if bond.GetBondType() == Chem.BondType.SINGLE:
                    next_atom = bond.GetOtherAtom(substituent_oxygen)
                    if next_atom.GetAtomicNum() == 6:
                        for next_bond in next_atom.GetBonds():
                            if next_bond.GetBondType() == Chem.BondType.DOUBLE:
                                other_atom = next_bond.GetOtherAtom(next_atom)
                                if other_atom.GetAtomicNum() == 6:
                                    vinyl_ether = True
                                    break

        substituent_type = "alkyl"  # default to alkyl if not ester or vinyl
        if ester:
            substituent_type = "acyl"
        elif vinyl_ether:
            substituent_type = "alk-1-enyl"

        return True, f"Glycerol with one {substituent_type} substituent"

    return False, "No valid substituent type found"