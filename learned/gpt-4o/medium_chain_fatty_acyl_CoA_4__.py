"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Updated CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)C[C@H](O)C(C)(C)COP(O)(=O)OCC1OC(n2cnc3nc[nH]c(N)c23)[C@@H](O)C1OP(O)(O)=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # Generalized thioester pattern for flexibility
    thioester_pattern = Chem.MolFromSmarts("C(=O)SC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for medium-chain length in the fatty acyl portion
    acyl_pattern = Chem.MolFromSmarts("C(=O)SC")
    matches = mol.GetSubstructMatches(acyl_pattern)
    medium_chain_found = False

    for match in matches:
        # Check chain length from the carbonyl adjacent to thioester
        carbon_count = 0
        atom_queue = [mol.GetAtomWithIdx(match[-1])]
        visited = set()
        
        while atom_queue:
            current_atom = atom_queue.pop(0)
            if current_atom.GetIdx() not in visited:
                visited.add(current_atom.GetIdx())
                if current_atom.GetSymbol() == 'C':
                    carbon_count += 1
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited and neighbor.GetSymbol() in ['C', 'H']:
                        atom_queue.append(neighbor)

        if 8 <= carbon_count <= 12:
            medium_chain_found = True
            break
    
    if not medium_chain_found:
        return False, "No medium-chain (8-12 carbons) fatty acyl portion found"

    # Check for the number of deprotonated oxygens
    deprotonated_oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if deprotonated_oxygen_count < 4:
        return False, f"Insufficient number of deprotonated oxygens: {deprotonated_oxygen_count}"

    return True, "Molecule is a medium-chain fatty acyl-CoA(4-)"