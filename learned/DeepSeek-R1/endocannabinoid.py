"""
Classifies: CHEBI:67197 endocannabinoid
"""
"""
Classifies: CHEBI:67194 endocannabinoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_endocannabinoid(smiles: str):
    """
    Determines if a molecule is an endocannabinoid based on its SMILES string.
    Endocannabinoids typically consist of a fatty acid derivative connected to a polar head group
    (ethanolamine or glycerol) via amide/ester bonds, with long carbon chains (C16+).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an endocannabinoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ethanolamide structure (fatty acid amide of ethanolamine)
    ethanolamide_pattern = Chem.MolFromSmarts("[CX3](=O)-[NX3]-[CX4]-[CX4]-[OH]")
    ethanolamide_matches = mol.GetSubstructMatches(ethanolamide_pattern)
    
    if ethanolamide_matches:
        # Check fatty acid chain length from amide carbonyl
        for match in ethanolamide_matches:
            carbonyl_atom = mol.GetAtomWithIdx(match[0])
            # Find attached carbon chain (non-amide/non-oxygen neighbor)
            chain_start = None
            for neighbor in carbonyl_atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(
                    carbonyl_atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                    continue
                if neighbor.GetIdx() == match[1]:  # Amide nitrogen
                    continue
                chain_start = neighbor
                break
            
            if chain_start and chain_start.GetAtomicNum() == 6:
                # Traverse chain to count carbons
                visited = set()
                stack = [(chain_start, 0)]  # (atom, depth)
                max_chain = 0
                while stack:
                    atom, depth = stack.pop()
                    if atom.GetAtomicNum() != 6 or atom in visited:
                        continue
                    visited.add(atom)
                    max_chain = max(max_chain, depth + 1)  # +1 to count current atom
                    for nbr in atom.GetNeighbors():
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                        if bond.GetBondType() == Chem.BondType.SINGLE:
                            stack.append((nbr, depth + 1))
                
                if max_chain >= 16:  # Minimum C16 chain
                    return True, "Ethanolamide structure with long fatty acid chain (C16+)"

    # Check for monoacylglycerol structure (glycerol + 1 ester)
    glycerol_ester_pattern = Chem.MolFromSmarts(
        "[CH2]-[CH](-O-C(=O)-[!O])-[CH2]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_ester_pattern)
    
    if glycerol_matches:
        # Verify glycerol hydroxyls and fatty acid chain length
        for match in glycerol_matches:
            # Check adjacent CH2 groups have hydroxyls
            left_ch2 = mol.GetAtomWithIdx(match[0])
            right_ch2 = mol.GetAtomWithIdx(match[5])
            has_left_oh = any(n.GetAtomicNum() == 8 for n in left_ch2.GetNeighbors())
            has_right_oh = any(n.GetAtomicNum() == 8 for n in right_ch2.GetNeighbors())
            
            if has_left_oh and has_right_oh:
                # Check ester chain length
                ester_carbonyl = mol.GetAtomWithIdx(match[3])
                chain_start = None
                for nbr in ester_carbonyl.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(
                        ester_carbonyl.GetIdx(), nbr.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                        continue
                    chain_start = nbr
                    break
                
                if chain_start and chain_start.GetAtomicNum() == 6:
                    # Traverse chain
                    visited = set()
                    stack = [(chain_start, 0)]
                    max_chain = 0
                    while stack:
                        atom, depth = stack.pop()
                        if atom.GetAtomicNum() != 6 or atom in visited:
                            continue
                        visited.add(atom)
                        max_chain = max(max_chain, depth + 1)
                        for nbr in atom.GetNeighbors():
                            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                            if bond.GetBondType() == Chem.BondType.SINGLE:
                                stack.append((nbr, depth + 1))
                    
                    if max_chain >= 16:  # Minimum C16 chain
                        return True, "Monoacylglycerol structure with long fatty acid chain (C16+)"

    # Additional checks could include dopamine conjugates or other head groups
    # but ethanolamides and monoacylglycerols cover major classes
    
    return False, "Does not match ethanolamide or monoacylglycerol patterns with C16+ chains"