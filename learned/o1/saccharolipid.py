"""
Classifies: CHEBI:166828 saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Find sugar rings (5 or 6-membered rings with exactly one oxygen)
    ssr = Chem.GetSymmSSSR(mol)
    sugar_rings = []
    for ring in ssr:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        ring_size = len(ring_atoms)
        if ring_size not in [5,6]:
            continue
        num_oxygen = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
        if num_oxygen != 1:
            continue
        # Check that all ring atoms are sp3 hybridized
        ring_ok = True
        for atom in ring_atoms:
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                ring_ok = False
                break
        if not ring_ok:
            continue
        sugar_rings.append(ring)
    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety found"

    # For each sugar ring, check for ester or amide linkages
    has_linkage = False
    for ring in sugar_rings:
        ring_atom_indices = list(ring)
        ring_atom_set = set(ring_atom_indices)
        for idx in ring_atom_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:
                continue
            # Look for oxygen or nitrogen connected to this carbon (potential linkage site)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring_atom_set:
                    continue
                if neighbor.GetAtomicNum() == 8:
                    # Check for ester linkage
                    # Neighbor oxygen should be connected to a carbonyl carbon
                    for nei_nei in neighbor.GetNeighbors():
                        if nei_nei.GetIdx() == atom.GetIdx():
                            continue
                        if nei_nei.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(neighbor.GetIdx(), nei_nei.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                            # Check if nei_nei has a double bond to oxygen (carbonyl group)
                            has_carbonyl = False
                            for nn_neighbor in nei_nei.GetNeighbors():
                                if nn_neighbor.GetIdx() == neighbor.GetIdx():
                                    continue
                                bond = mol.GetBondBetweenAtoms(nei_nei.GetIdx(), nn_neighbor.GetIdx())
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nn_neighbor.GetAtomicNum() == 8:
                                    has_carbonyl = True
                                    break
                            if has_carbonyl:
                                has_linkage = True
                                break
                elif neighbor.GetAtomicNum() == 7:
                    # Check for amide linkage
                    # Neighbor nitrogen should be connected to a carbonyl carbon
                    for nei_nei in neighbor.GetNeighbors():
                        if nei_nei.GetIdx() == atom.GetIdx():
                            continue
                        if nei_nei.GetAtomicNum() == 6 and mol.GetBondBetweenAtoms(neighbor.GetIdx(), nei_nei.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                            # Check if nei_nei has a double bond to oxygen (carbonyl group)
                            has_carbonyl = False
                            for nn_neighbor in nei_nei.GetNeighbors():
                                if nn_neighbor.GetIdx() == neighbor.GetIdx():
                                    continue
                                bond = mol.GetBondBetweenAtoms(nei_nei.GetIdx(), nn_neighbor.GetIdx())
                                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nn_neighbor.GetAtomicNum() == 8:
                                    has_carbonyl = True
                                    break
                            if has_carbonyl:
                                has_linkage = True
                                break
                if has_linkage:
                    break
            if has_linkage:
                break
        if has_linkage:
            break
    if not has_linkage:
        return False, "No ester or amide linkage found between lipid and sugar moieties"

    return True, "Contains both carbohydrate and lipid moieties linked via ester or amide bonds"