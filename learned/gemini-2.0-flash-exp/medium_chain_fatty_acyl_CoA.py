"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA consists of Coenzyme A linked to a fatty acid with 6-12 carbons via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for CoA core (Adenine-Ribose-Phosphate) - More specific pattern
    core_pattern = Chem.MolFromSmarts("n1cnc2c(N)ncnc12[C@H]1[C@@H]([C@@H]([C@H](O[P]([OX1])(=[OX1])O)O)O)[C@H]1COP([OX1])(=[OX1])OP([OX1])(=[OX1])O")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "CoA core (Adenine-Ribose-Phosphate) not found"

    # 2. Check for pantetheine-like moiety
    pantetheine_pattern = Chem.MolFromSmarts("NCCSC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(pantetheine_pattern):
        return False, "Pantetheine-like chain not found"
    
    # 3. Check for thioester linkage (-C(=O)-S-)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester groups, need exactly 1"
    
    # 4. Analyze fatty acid chain length
    for match in thioester_matches:
        carbonyl_c_index = match[0]
        thio_s_index = match[1]
        
        carbonyl_c_atom = mol.GetAtomWithIdx(carbonyl_c_index)
        
        # Find the carbon directly attached to the carbonyl carbon of the thioester (start of fatty acid chain)
        connected_carbon_indices = [neighbor.GetIdx() for neighbor in carbonyl_c_atom.GetNeighbors() if neighbor.GetIdx() != thio_s_index and neighbor.GetAtomicNum() == 6]

        if len(connected_carbon_indices) != 1:
            return False, "Could not identify the connecting carbon to the carbonyl in thioester."

        start_carbon_index = connected_carbon_indices[0]

        # Traverse the fatty acid chain and count heavy atoms
        chain_length = 0
        visited_atoms = set()
        current_atom_index = start_carbon_index
        prev_atom_index = carbonyl_c_index
        
        while True:
            visited_atoms.add(current_atom_index)
            current_atom = mol.GetAtomWithIdx(current_atom_index)
            is_terminal = True
            next_atom_index = None
            for neighbor in current_atom.GetNeighbors():
                neighbor_index = neighbor.GetIdx()
                if neighbor.GetAtomicNum() != 1 and neighbor_index != prev_atom_index and neighbor_index not in visited_atoms:
                    next_atom_index = neighbor_index
                    is_terminal = False
                    break
            if is_terminal:
                 break
            if not current_atom.GetAtomicNum()==6 and current_atom.GetAtomicNum() != 1:
                chain_length += 1 #Count only heavy atoms
            
            prev_atom_index = current_atom_index
            current_atom_index = next_atom_index

        # Count last atom if heavy atom
        last_atom = mol.GetAtomWithIdx(current_atom_index)
        if last_atom.GetAtomicNum() != 1:
          chain_length +=1
            
        if not 6 <= chain_length <= 12:
            return False, f"Fatty acid chain has {chain_length} heavy atoms, should be 6-12."


    return True, "Medium-chain fatty acyl-CoA identified"