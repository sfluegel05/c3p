"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: branched-chain fatty acyl-CoA
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA is a fatty acyl-CoA resulting from the condensation 
    of the thiol group of coenzyme A with the carboxy group of any branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    from rdkit import Chem
    from collections import deque

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thioester linkage pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if not thioester_matches:
        return False, "No thioester linkage found"

    # Check if the sulfur is connected to the CoA moiety
    # Define a simplified CoA pattern focusing on the pantetheine moiety
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Sulfur is not connected to CoA moiety"

    # Identify the fatty acyl chain attached via thioester bond
    # Find the carbonyl carbon and sulfur atom of the thioester linkage
    for match in thioester_matches:
        carbonyl_c_idx = match[0]
        sulfur_idx = match[2]
        # Get the fatty acyl chain starting from carbonyl carbon
        fatty_acyl_atoms = get_acyl_chain(mol, carbonyl_c_idx, sulfur_idx)
        if fatty_acyl_atoms:
            break
    if not fatty_acyl_atoms:
        return False, "Could not identify fatty acyl chain"

    # Check for branching in the fatty acyl chain
    branching_points = 0
    for atom_idx in fatty_acyl_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon atom
            neighbor_carbons_in_chain = [nbr.GetIdx() for nbr in atom.GetNeighbors() 
                                         if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in fatty_acyl_atoms]
            if len(neighbor_carbons_in_chain) > 2:
                branching_points += 1
    if branching_points == 0:
        return False, "Fatty acyl chain is not branched"
    else:
        return True, f"Fatty acyl chain is branched with {branching_points} branching point(s)"

def get_acyl_chain(mol, start_atom_idx, exclude_atom_idx):
    """
    Traverses the molecule to get the fatty acyl chain starting from the carbonyl carbon.

    Args:
        mol (rdkit.Chem.Mol): Molecule object
        start_atom_idx (int): Index of the carbonyl carbon atom
        exclude_atom_idx (int): Index of the sulfur atom to exclude from traversal

    Returns:
        List[int]: List of atom indices in the fatty acyl chain
    """
    from collections import deque

    visited = set()
    queue = deque([start_atom_idx])
    acyl_chain = []

    while queue:
        atom_idx = queue.popleft()
        if atom_idx in visited:
            continue
        visited.add(atom_idx)
        acyl_chain.append(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            # Exclude sulfur atom and atoms connected to sulfur
            if nbr_idx == exclude_atom_idx:
                continue
            if neighbor.GetAtomicNum() not in [1, 6, 8]:  # Only traverse C, H, O atoms
                continue
            queue.append(nbr_idx)
    return acyl_chain

__metadata__ = {  
    'chemical_class': {  
        'name': 'branched-chain fatty acyl-CoA',
        'definition': 'A fatty acyl-CoA that results from the formal condensation of the thiol group of coenzyme A with the carboxy group of any branched-chain fatty acid.'},
    'message': None,
    'success': True,
    'error': '',
    'stdout': None
}