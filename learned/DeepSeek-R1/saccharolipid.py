"""
Classifies: CHEBI:166828 saccharolipid
"""
"""
Classifies: CHEBI: saccharolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_saccharolipid(smiles: str):
    """
    Determines if a molecule is a saccharolipid based on its SMILES string.
    Saccharolipids contain both carbohydrate (sugar) and lipid components covalently linked.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a saccharolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbohydrate component: cyclic structure with multiple hydroxyl groups
    # Look for a 5 or 6-membered ring with at least 3 oxygen atoms (as hydroxyl or ring oxygen)
    carb_pattern = Chem.MolFromSmarts('[C;R5,R6][OH]')
    carb_matches = mol.GetSubstructMatches(carb_pattern)
    if len(carb_matches) < 3:
        return False, "Insufficient carbohydrate features (cyclic structure with hydroxyls)"

    # Check for lipid component: long aliphatic chain (>=8 carbons) or ester/amide linkage
    # Look for ester (O-C=O) or amide (N-C=O) groups connected to long chains
    lipid_found = False
    ester_amide_pattern = Chem.MolFromSmarts('[O,N][C]=O')
    ester_amide_matches = mol.GetSubstructMatches(ester_amide_pattern)
    
    for match in ester_amide_matches:
        # Get the oxygen/nitrogen atom in the ester/amide
        o_n_atom = mol.GetAtomWithIdx(match[0])
        # Check neighboring carbons for long chains
        for neighbor in o_n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                # Traverse the chain starting from this carbon
                chain_length = 0
                current = neighbor
                visited = set()
                stack = [(current, 0)]  # (atom, current_length)
                max_length = 0
                while stack:
                    atom, length = stack.pop()
                    if atom.GetIdx() in visited:
                        continue
                    visited.add(atom.GetIdx())
                    if length > max_length:
                        max_length = length
                    # Follow non-ring, non-carbonyl carbons
                    for next_atom in atom.GetNeighbors():
                        if (next_atom.GetAtomicNum() == 6 and 
                            not next_atom.IsInRing() and 
                            next_atom.GetBondToAtom(atom).GetBondType() == Chem.BondType.SINGLE):
                            stack.append((next_atom, length + 1))
                if max_length >= 7:  # 8 carbons including starting point
                    lipid_found = True
                    break
        if lipid_found:
            break

    if not lipid_found:
        return False, "No lipid component (long chain with ester/amide) detected"

    return True, "Contains carbohydrate and lipid components linked via ester/amide"