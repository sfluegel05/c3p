"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolFromSmarts
from rdkit.Chem.rdmolops import GetAdjacencyMatrix

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside.
    A pyrimidine deoxyribonucleoside consists of a deoxyribose sugar connected to a pyrimidine base via a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Check for deoxyribose sugar: five-membered ring with one oxygen and a missing hydroxyl at C2'
    sugar_pattern = MolFromSmarts("[O;R]1[C;R][C;R][C;R][C;R]1")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No five-membered sugar ring with oxygen found"
    
    deoxyribose_found = False
    for sugar_match in sugar_matches:
        # Get the oxygen atom in the sugar ring
        oxygen_idx = [atom for atom in sugar_match if mol.GetAtomWithIdx(atom).GetAtomicNum() == 8][0]
        oxygen = mol.GetAtomWithIdx(oxygen_idx)
        
        # Find adjacent carbon in the ring (C1')
        neighbors = [n for n in oxygen.GetNeighbors() if n.GetIdx() in sugar_match]
        if not neighbors:
            continue
        c1 = neighbors[0]
        
        # Check C2' (next carbon in the ring) has no hydroxyl
        c1_neighbors = [n for n in c1.GetNeighbors() if n.GetIdx() in sugar_match]
        if len(c1_neighbors) < 2:
            continue
        next_c = [n for n in c1_neighbors if n.GetIdx() != oxygen_idx][0]
        if next_c.GetTotalNumHs() < 2 or any(bond.GetBondType() == Chem.BondType.SINGLE and 
                                            bond.GetOtherAtom(next_c).GetAtomicNum() == 8 for bond in next_c.GetBonds()):
            continue  # C2' has a hydroxyl or other oxygen substituent
        
        # Check for hydroxymethyl group (C5')
        c5_candidates = [atom for atom in sugar_match if mol.GetAtomWithIdx(atom).GetAtomicNum() == 6 and
                         any(n.GetAtomicNum() == 8 and n.GetTotalNumHs() > 0 for n in mol.GetAtomWithIdx(atom).GetNeighbors())]
        if not c5_candidates:
            continue
        
        deoxyribose_found = True
        break
    
    if not deoxyribose_found:
        return False, "No deoxyribose sugar found"
    
    # Check for pyrimidine base (six-membered ring with two nitrogens)
    pyrimidine_pattern = MolFromSmarts("n1cnccc1")  # Basic pyrimidine skeleton
    base_matches = mol.GetSubstructMatches(pyrimidine_pattern)
    if not base_matches:
        return False, "No pyrimidine base detected"
    
    # Check glycosidic bond between sugar and base
    glycosidic_bond = False
    for base_match in base_matches:
        for atom_idx in base_match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() == 7:  # Nitrogen in base
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in sugar_match:
                        glycosidic_bond = True
                        break
                if glycosidic_bond:
                    break
        if glycosidic_bond:
            break
    
    if not glycosidic_bond:
        return False, "No glycosidic bond between sugar and base"
    
    # Check for absence of phosphate groups (nucleoside, not nucleotide)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Phosphate group present"
    
    return True, "Deoxyribose sugar with pyrimidine base connected via glycosidic bond"