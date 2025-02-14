"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Draw


def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    Fatty acids are characterized by a long aliphatic chain and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
         return False, "Must have at least one carboxylic acid group"
         
    # 2. Check for common atoms and a main chain composed of C,H and possibly O.
    allowed_atoms = [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]  # H, C, N, O, F, P, S, Cl, Br, I
    main_chain_atoms = [6,1,8]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Non-allowed atoms present"

    # 3. Check number of carbons - should have at least 4 in the main chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 4:
         return False, "Chain too short to be a fatty acid"
    
    # 4. Analyze the longest chain (as a proxy for main chain)
    longest_chain_atoms = []
    max_chain_len = 0;
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom().GetIdx()
        a2 = bond.GetEndAtom().GetIdx()
        
        if bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6:
             
                visited = []
                
                def find_chain(current_atom_index, chain, visited):
                    visited.append(current_atom_index)
                    max_len = len(chain);
                    max_chain = chain.copy()

                    neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(current_atom_index).GetNeighbors()]
                    
                    for neighbor_idx in neighbors:
                        neighbor = mol.GetAtomWithIdx(neighbor_idx)
                        if neighbor.GetAtomicNum() == 6 and neighbor_idx not in visited:
                            
                             new_chain, new_max_len = find_chain(neighbor_idx, chain+[neighbor_idx], visited.copy());
                             if new_max_len > max_len:
                                max_len = new_max_len;
                                max_chain = new_chain
                    return max_chain, max_len;

                chain1, len1 = find_chain(a1, [a1], visited.copy())
                if len1 > max_chain_len:
                     max_chain_len = len1;
                     longest_chain_atoms = chain1

                chain2, len2 = find_chain(a2, [a2], visited.copy())
                if len2 > max_chain_len:
                     max_chain_len = len2;
                     longest_chain_atoms = chain2

    if max_chain_len < 4:
         return False, "Chain too short to be a fatty acid"
    
    # 5. Check for branches (max 4 non-hydrogen branch points on the main chain) - use all other atoms as branch points
    non_main_chain_atoms = 0;
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1 and atom.GetIdx() not in longest_chain_atoms:
             non_main_chain_atoms +=1;
    if non_main_chain_atoms > 4:
         return False, "Too many branches for a typical fatty acid"

    # 6. Check for rings that are *part of the main chain*. If rings are found, the atoms must be part of the chain.
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
         for ring in ring_info.AtomRings():
              is_chain_ring = True;
              for atom_index in ring:
                   if atom_index not in longest_chain_atoms:
                        is_chain_ring = False;
                        break;
              if not is_chain_ring:
                    return False, "Rings cannot be outside the main carbon chain."

    return True, "Meets criteria for a fatty acid"