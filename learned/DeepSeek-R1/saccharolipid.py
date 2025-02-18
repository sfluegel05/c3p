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

    # Improved carbohydrate detection: look for sugar-like rings (5/6 membered with multiple oxygens)
    # Pattern matches rings with >=2 oxygen atoms (either in ring or as hydroxyls)
    carb_pattern = Chem.MolFromSmarts('[O;r5,r6]~[C;r5,r6]~[O;r5,r6]')
    if not mol.HasSubstructMatch(carb_pattern):
        carb_pattern2 = Chem.MolFromSmarts('[C;r5,r6][OH]')
        carb_matches = mol.GetSubstructMatches(carb_pattern2)
        if len(carb_matches) < 2:
            return False, "No carbohydrate moiety detected (cyclic structure with multiple oxygens)"

    # Lipid component check: long aliphatic chain (>=8 carbons) connected via ester/amide
    lipid_found = False
    ester_amide_pattern = Chem.MolFromSmarts('[O,N][C]=O')
    ester_amide_matches = mol.GetSubstructMatches(ester_amide_pattern)
    
    for match in ester_amide_matches:
        # Start from carbonyl carbon (index 1 in match [O,N]-C=O)
        carbonyl_carbon = mol.GetAtomWithIdx(match[1])
        
        # Traverse two directions from carbonyl carbon
        for neighbor in carbonyl_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() != 6 or neighbor.GetIdx() == match[0]:
                continue  # Skip non-carbon atoms and back to O/N
                
            # Track longest chain in this direction
            chain_length = 0
            visited = set()
            stack = [(neighbor, 0)]
            
            while stack:
                atom, depth = stack.pop()
                if atom.GetIdx() in visited:
                    continue
                visited.add(atom.GetIdx())
                
                # Only follow aliphatic carbons not in rings
                if (atom.GetAtomicNum() == 6 and
                    not atom.IsInRing() and
                    not any(bond.GetBondType() == Chem.BondType.AROMATIC 
                            for bond in atom.GetBonds())):
                    
                    chain_length = max(chain_length, depth + 1)
                    # Continue traversing non-ring, non-carbonyl carbons
                    for next_atom in atom.GetNeighbors():
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), next_atom.GetIdx())
                        if (bond.GetBondType() == Chem.BondType.SINGLE and
                            next_atom.GetAtomicNum() == 6 and
                            not next_atom.IsInRing()):
                            stack.append((next_atom, depth + 1))
            
            if chain_length >= 7:  # 8 carbons including starting point
                lipid_found = True
                break
        if lipid_found:
            break

    if not lipid_found:
        return False, "No lipid component (long chain >=8 carbons connected via ester/amide)"

    return True, "Contains carbohydrate and lipid components linked via ester/amide"