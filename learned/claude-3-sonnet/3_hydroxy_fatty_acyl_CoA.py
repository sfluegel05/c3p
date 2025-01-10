"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
"""
Classifies: 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for adenine core
    adenine = Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(adenine):
        return False, "No adenine moiety found"

    # Look for ribose with phosphate
    ribose_phosphate = Chem.MolFromSmarts("OCC1OC(n2cnc3c2ncnc3)C(O)C1OP(O)(=O)O")
    if not mol.HasSubstructMatch(ribose_phosphate):
        return False, "Missing ribose-phosphate structure"

    # Check for pantetheine arm with thioester
    pantetheine = Chem.MolFromSmarts("CC(C)(COP(O)(=O)OP(O)(=O)O)[C@H](O)C(=O)NCCC(=O)NCCS")
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Missing or incomplete pantetheine arm"

    # Check for thioester linkage
    thioester = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester):
        return False, "No thioester linkage found"

    # Check for 3-hydroxy pattern near thioester
    # Match both R and S configurations, including unspecified stereochemistry
    hydroxy_patterns = [
        "[CH2X4]-[CH1X4]([OX2H])-[CH2X4]-C(=O)S", # Generic
        "[CH2X4]-[C@H]([OX2H])-[CH2X4]-C(=O)S",   # S config
        "[CH2X4]-[C@@H]([OX2H])-[CH2X4]-C(=O)S"   # R config
    ]
    
    has_hydroxy = False
    for pattern in hydroxy_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            has_hydroxy = True
            break
            
    if not has_hydroxy:
        return False, "No 3-hydroxy group found adjacent to thioester"

    # Count phosphate groups (including various oxidation states)
    phosphate_pattern = Chem.MolFromSmarts("P(~O)(~O)(~O)~O")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_count < 3:
        return False, f"Insufficient phosphate groups (found {phosphate_count}, need at least 3)"

    # Verify presence of fatty acid chain by counting carbons between OH and thioester
    # First find the thioester carbon
    thioester_matches = mol.GetSubstructMatches(thioester)
    if not thioester_matches:
        return False, "Could not locate thioester group"
        
    # Find hydroxy group connected to carbon chain
    hydroxy_carbon = None
    for pattern in hydroxy_patterns:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        if matches:
            hydroxy_carbon = matches[0][1]  # Index of carbon bearing OH
            break
    
    if hydroxy_carbon is None:
        return False, "Could not locate 3-hydroxy carbon"

    # Check if we have a reasonable fatty acid chain
    # Count carbons in main chain using atomic numbers and connectivity
    carbon_count = 0
    visited = set()
    def count_chain_carbons(atom_idx, prev_idx=None):
        if atom_idx in visited:
            return
        visited.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:  # Carbon
            nonlocal carbon_count
            carbon_count += 1
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() != prev_idx and neighbor.GetAtomicNum() == 6:
                    count_chain_carbons(neighbor.GetIdx(), atom_idx)

    count_chain_carbons(hydroxy_carbon)
    
    if carbon_count < 4:
        return False, f"Fatty acid chain too short (found approximately {carbon_count} carbons)"

    return True, "Contains CoA moiety with 3-hydroxy fatty acid chain attached via thioester linkage"