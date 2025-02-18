"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:XXXXX medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if matches criteria, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Match CoA core structure with stereochemistry
    coa_core = Chem.MolFromSmarts(
        "[NH]C(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12"
    )
    if not mol.HasSubstructMatch(coa_core):
        return False, "Missing CoA core structure"

    # Find thioester group (S-C(=O)-R)
    thioester = Chem.MolFromSmarts("[SX2][CX3](=O)")
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester bond found"

    sulfur_idx = thioester_matches[0][0]
    carbonyl_atom = mol.GetAtomWithIdx(thioester_matches[0][1])

    def find_aliphatic_chain(atom, visited=None):
        if visited is None:
            visited = set()
        if atom.GetIdx() in visited or atom.GetAtomicNum() != 6:
            return 0
        if atom.IsInRing():
            return 0
        visited.add(atom.GetIdx())
        
        max_length = 0
        for neighbor in atom.GetNeighbors():
            # Prevent backtracking to CoA through sulfur
            if neighbor.GetIdx() == sulfur_idx:
                continue
            # Only follow aliphatic carbon chains
            if neighbor.GetAtomicNum() != 6 or neighbor.IsInRing():
                continue
            
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
            # Allow single or double bonds between carbons
            if bond.GetBondType() not in (Chem.BondType.SINGLE, Chem.BondType.DOUBLE):
                continue
            
            length = find_aliphatic_chain(neighbor, visited.copy())
            max_length = max(max_length, length)
            
        return 1 + max_length

    chain_length = find_aliphatic_chain(carbonyl_atom)
    if not (6 <= chain_length <= 12):
        return False, f"Aliphatic acyl chain length {chain_length} not in medium range (6-12)"

    # Verify total charge is -4 (3 phosphate groups deprotonated + thioester)
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != -4:
        return False, f"Total charge {total_charge} ≠ -4"

    return True, "Medium-chain fatty acyl-CoA(4-) with correct structure and charge"