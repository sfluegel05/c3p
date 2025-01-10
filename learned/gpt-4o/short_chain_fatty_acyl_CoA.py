"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

    A short-chain fatty acyl-CoA is a fatty acyl-CoA formed by the condensation
    of the thiol group of coenzyme A with a short-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Redefine Coenzyme A with a more comprehensive and flexible pattern
    coenzymeA_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)C[S]CCNC(=O)C[C@H](O)COP(O)(=O)OP(O)(O)=OOC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c(ncnc12)N")
    if not mol.HasSubstructMatch(coenzymeA_pattern):
        return False, "No Coenzyme A moiety found"
    
    # Detect thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)S[C]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Ensure short-chain fatty acid component (2-6 carbons)
    # More flexible carbon counting logic
    for match in mol.GetSubstructMatches(thioester_pattern):
        sulfur_connected_c = [atom.GetIdx() for atom in mol.GetAtomWithIdx(match[1]).GetNeighbors() if atom.GetSymbol() == 'C']
        if not sulfur_connected_c:
            continue
        start_atom_idx = sulfur_connected_c[0]
        
        visited = set()
        carbon_count = 0
        to_visit = [start_atom_idx]
        while to_visit:
            current_atom_idx = to_visit.pop()
            if current_atom_idx not in visited:
                current_atom = mol.GetAtomWithIdx(current_atom_idx)
                if current_atom.GetSymbol() == 'C':
                    carbon_count += 1
                    visited.add(current_atom_idx)
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetIdx() not in visited and neighbor.GetSymbol() == 'C':
                            to_visit.append(neighbor.GetIdx())
        
        if 2 <= carbon_count <= 6:
            return True, "Contains a CoA moiety joined by thioester bond to a short-chain fatty acid"
    
    return False, "No proper short-chain fatty acyl component found"