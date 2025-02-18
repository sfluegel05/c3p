"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:XXXXX medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors

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

    # Verify CoA core structure components
    coa_core = Chem.MolFromSmarts(
        "[NH]C(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC([C@H]([C@@H]1O)O)n1cnc2c(N)ncnc12"
    )
    if not mol.HasSubstructMatch(coa_core):
        return False, "Missing CoA core structure"

    # Check for thioester-linked acyl group (S-C(=O)-R)
    thioester = Chem.MolFromSmarts("[SX2][CX3](=O)")
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "No thioester bond found"

    # Get the sulfur and carbonyl atoms from the first thioester match
    sulfur_idx = thioester_matches[0][0]
    carbonyl_idx = thioester_matches[0][1]
    sulfur = mol.GetAtomWithIdx(sulfur_idx)
    carbonyl = mol.GetAtomWithIdx(carbonyl_idx)

    # Traverse the acyl chain starting from the carbonyl carbon
    chain_carbons = set()
    stack = [carbonyl]
    visited = set()

    while stack:
        current = stack.pop()
        if current.GetIdx() in visited:
            continue
        visited.add(current.GetIdx())
        if current.GetAtomicNum() != 6:  # Only consider carbon atoms
            continue
        chain_carbons.add(current.GetIdx())
        # Traverse all neighbors except sulfur and carbonyl oxygen
        for neighbor in current.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            # Skip sulfur and carbonyl oxygen
            if neighbor_idx == sulfur.GetIdx():
                continue
            # Skip oxygen double-bonded to carbonyl (part of thioester)
            if neighbor.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(current.GetIdx(), neighbor_idx).GetBondType() == Chem.BondType.DOUBLE:
                continue
            # Skip hydrogens
            if neighbor.GetAtomicNum() == 1:
                continue
            if neighbor_idx not in visited:
                stack.append(neighbor)

    chain_length = len(chain_carbons)
    if not (6 <= chain_length <= 12):
        return False, f"Acyl chain length {chain_length} not in medium range (6-12)"

    # Check for four deprotonated phosphate oxygens
    phosphate_neg = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 15:  # Phosphorus
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 8 and neighbor.GetFormalCharge() == -1:
                    phosphate_neg += 1
    if phosphate_neg != 4:
        return False, f"Found {phosphate_neg} deprotonated phosphate oxygens, need 4"

    return True, "Medium-chain fatty acyl-CoA(4-) with correct structure and charge"