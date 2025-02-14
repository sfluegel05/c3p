"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    A saccharolipid is defined as a lipid that contains a carbohydrate moiety linked via ester or amide bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    sugar_pattern = Chem.MolFromSmarts("[#6R1]-[@#6R1]-[@#6R1]-[@#6R1]-[@#6R1]-[@#8R1]")  # Pyranose ring
    furanose_pattern = Chem.MolFromSmarts("[#6R1]-[@#6R1]-[@#6R1]-[@#6R1]-[@#8R1]")      # Furanose ring
    long_chain_pattern = Chem.MolFromSmarts("C[C@@H](CC)CC[C@@H](C)CCCCCCCC")            # Chain of at least 8 carbons
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    amide_pattern = Chem.MolFromSmarts("C(=O)N")

    # Find sugar rings
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    sugar_rings = sugar_matches + furanose_matches

    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety found"

    # Find long aliphatic chains
    chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if not chain_matches:
        return False, "No lipid (long aliphatic chain) moiety found"

    # Find ester or amide linkages
    ester_bonds = []
    amide_bonds = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
            continue
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        begin_idx = begin_atom.GetIdx()
        end_idx = end_atom.GetIdx()
        atom_pair = {begin_idx, end_idx}

        if mol.HasSubstructMatch(Chem.MolFromSmarts(f"[#{begin_atom.GetAtomicNum()}]({ester_pattern})")) and \
           mol.HasSubstructMatch(Chem.MolFromSmarts(f"[#{end_atom.GetAtomicNum()}]({ester_pattern})")):
            ester_bonds.append(bond)
        if mol.HasSubstructMatch(Chem.MolFromSmarts(f"[#{begin_atom.GetAtomicNum()}]({amide_pattern})")) and \
           mol.HasSubstructMatch(Chem.MolFromSmarts(f"[#{end_atom.GetAtomicNum()}]({amide_pattern})")):
            amide_bonds.append(bond)

    if not ester_bonds and not amide_bonds:
        return False, "No ester or amide linkage found between lipid and sugar moieties"

    # Check if ester/amide linkages connect sugar rings and lipid chains
    has_linkage = False
    for bond in ester_bonds + amide_bonds:
        begin_atom = bond.GetBeginAtom()
        end_atom = bond.GetEndAtom()
        begin_idx = begin_atom.GetIdx()
        end_idx = end_atom.GetIdx()

        # Check if one atom is in sugar ring and the other is in lipid chain
        in_sugar = False
        in_chain = False
        for ring in sugar_rings:
            if begin_idx in ring or end_idx in ring:
                in_sugar = True
                break
        for chain in chain_matches:
            if begin_idx in chain or end_idx in chain:
                in_chain = True
                break
        if in_sugar and in_chain:
            has_linkage = True
            break

    if not has_linkage:
        return False, "Ester or amide bonds do not link sugar and lipid moieties"

    return True, "Contains carbohydrate and lipid moieties linked via ester or amide bonds"