"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:XXXXX long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    A long-chain fatty acyl-CoA is a thioester formed between coenzyme A and a long-chain fatty acid (C13-C22).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA structural features
    # 1. Pantetheine fragment (S-C-C-N-C(=O))
    pantetheine_pattern = Chem.MolFromSmarts('[S]-[C]-[C]-[N]-[C]=[O]')
    pantetheine_matches = mol.GetSubstructMatches(pantetheine_pattern)
    if not pantetheine_matches:
        return False, "Pantetheine fragment not found"
    
    # 2. At least two phosphate groups (indicative of CoA)
    phosphate_pattern = Chem.MolFromSmarts('[O]P(=O)([O])[O]')
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 2:
        return False, "Insufficient phosphate groups for CoA"
    
    # Find thioester groups (S-C(=O))
    thioester_pattern = Chem.MolFromSmarts('[S]-[C]=[O]')
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester group found"
    
    # Get all pantetheine sulfur atoms
    pantetheine_sulfur = {match[0] for match in pantetheine_matches}
    
    # Check each thioester group connected to pantetheine
    for match in thioester_matches:
        s_idx, c_idx = match[0], match[1]
        if s_idx not in pantetheine_sulfur:
            continue  # thioester not part of CoA
        
        # Get R group atom (connected to carbonyl carbon, not S)
        carbonyl_carbon = mol.GetAtomWithIdx(c_idx)
        r_group_atom = None
        for neighbor in carbonyl_carbon.GetNeighbors():
            if neighbor.GetIdx() == s_idx:
                continue
            if neighbor.GetAtomicNum() == 8:  # skip oxygen in carbonyl
                continue
            r_group_atom = neighbor
            break
        if not r_group_atom:
            continue
        
        # Traverse R group to count carbons and check validity
        visited = set()
        stack = [r_group_atom]
        carbon_count = 0
        valid_chain = True
        
        while stack and valid_chain:
            atom = stack.pop()
            if atom.GetIdx() in visited:
                continue
            visited.add(atom.GetIdx())
            
            # Check atom type
            atomic_num = atom.GetAtomicNum()
            if atomic_num not in [6, 8]:
                valid_chain = False
                break
            
            # Oxygen must be in carbonyl or hydroxyl
            if atomic_num == 8:
                oxygen_bonds = atom.GetBonds()
                valid_oxygen = False
                for bond in oxygen_bonds:
                    other_atom = bond.GetOtherAtom(atom)
                    if bond.GetBondType() == Chem.BondType.DOUBLE and other_atom.GetAtomicNum() == 6:
                        valid_oxygen = True  # carbonyl oxygen
                    elif bond.GetBondType() == Chem.BondType.SINGLE and other_atom.GetAtomicNum() == 6:
                        valid_oxygen = True  # hydroxyl oxygen
                if not valid_oxygen:
                    valid_chain = False
                    break
            
            # Check for aromatic rings in the chain
            if atom.GetIsAromatic():
                valid_chain = False
                break
            
            # Count carbons
            if atomic_num == 6:
                carbon_count += 1
            
            # Continue traversal
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() == c_idx or neighbor.GetAtomicNum() == 16:  # avoid backtracking
                    continue
                stack.append(neighbor)
        
        if not valid_chain:
            continue
        
        # Check chain length (C13-C22 corresponds to R group C12-C21)
        if 12 <= carbon_count <= 21:
            return True, f"Long-chain acyl group ({carbon_count + 1} carbons) attached via thioester to CoA"
    
    return False, "Acyl chain length not in C13-C22 range or invalid structure"