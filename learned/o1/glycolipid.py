"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: Glycolipid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as any molecule consisting of a glycosidic linkage
    between a carbohydrate moiety (usually a mono-, di-, or trisaccharide)
    and a lipid moiety (such as fatty acids, sphingolipids, or glycerolipids).
    In some cases, the glycerol backbone may be absent, and the sugar part may be acylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carbohydrate moiety
    has_sugar = False
    sugar_atoms = set()
    # Define general patterns for pyranose and furanose rings
    pyranose_pattern = Chem.MolFromSmarts("C1C[C@H]([O])[C@@H](O)[C@@H](O)O1")
    furanose_pattern = Chem.MolFromSmarts("C1C[C@H](O)[C@@H](O)O1")
    # Search for rings matching the sugar patterns
    matches = mol.GetSubstructMatches(pyranose_pattern) + mol.GetSubstructMatches(furanose_pattern)
    if matches:
        has_sugar = True
        for match in matches:
            sugar_atoms.update(match)
    else:
        # Alternatively, identify rings with oxygen and multiple hydroxyl groups
        ri = mol.GetRingInfo()
        for ring in ri.AtomRings():
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            if len(ring_atoms) in (5, 6):  # Five or six-membered rings
                o_in_ring = any(atom.GetSymbol() == 'O' for atom in ring_atoms)
                if o_in_ring:
                    # Count number of hydroxyl groups attached to ring carbons
                    num_oh = 0
                    for atom in ring_atoms:
                        if atom.GetSymbol() == 'C':
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetSymbol() == 'O' and neighbor.GetIdx() not in ring:
                                    # Check if oxygen is a hydroxyl group
                                    if neighbor.GetTotalDegree() == 1:
                                        num_oh += 1
                    if num_oh >= 2:
                        has_sugar = True
                        for atom in ring_atoms:
                            sugar_atoms.add(atom.GetIdx())
                        break
    if not has_sugar:
        return False, "No carbohydrate moiety found"

    # Identify lipid moiety
    has_lipid = False
    lipid_atoms = set()
    # Define a pattern for long aliphatic chains (length >=8)
    aliphatic_chain_pattern = Chem.MolFromSmarts("[C&R0][C&R0][C&R0][C&R0][C&R0][C&R0][C&R0][C&R0]")
    chain_matches = mol.GetSubstructMatches(aliphatic_chain_pattern)
    if chain_matches:
        has_lipid = True
        for match in chain_matches:
            lipid_atoms.update(match)
    else:
        # Also check for sphingoid bases
        sphingoid_pattern = Chem.MolFromSmarts("NC[C@H](O)CC=C[C@H](O)CC[C@H](O)CC")
        sphingoid_matches = mol.GetSubstructMatches(sphingoid_pattern)
        if sphingoid_matches:
            has_lipid = True
            for match in sphingoid_matches:
                lipid_atoms.update(match)
    if not has_lipid:
        return False, "No lipid moiety found"

    # Check for linkage between carbohydrate and lipid moieties
    linked = False
    for s_atom_idx in sugar_atoms:
        s_atom = mol.GetAtomWithIdx(s_atom_idx)
        for neighbor in s_atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx in lipid_atoms:
                linked = True
                break
        if linked:
            break

    if linked:
        return True, "Contains carbohydrate and lipid moieties linked together"
    else:
        # Also check for acylated sugars (sugar part acylated by fatty acids)
        for s_atom_idx in sugar_atoms:
            s_atom = mol.GetAtomWithIdx(s_atom_idx)
            if s_atom.GetSymbol() == 'O':
                for neighbor in s_atom.GetNeighbors():
                    if neighbor.GetSymbol() == 'C':
                        # Check if this carbon is part of a lipid chain
                        path = Chem.rdmolops.GetShortestPath(mol, neighbor.GetIdx(), list(lipid_atoms))
                        if path:
                            linked = True
                            break
                if linked:
                    break
        if linked:
            return True, "Contains carbohydrate moiety acylated by fatty acid"
        return False, "No linkage between carbohydrate and lipid moieties found"