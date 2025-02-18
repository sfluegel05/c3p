"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
"""
Classifies: 2-acyl-1-alkyl-sn-glycero-3-phosphocholine (CHEBI:XXXXX)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine.
    Criteria:
    - sn-glycerol backbone with correct stereochemistry
    - Phosphocholine group at sn-3
    - Ether-linked alkyl chain at sn-1
    - Ester-linked acyl chain at sn-2
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Core structure with chiral center, phosphocholine, ether, and ester
    core_pattern = Chem.MolFromSmarts("""
        [C@H](COP(=O)([O-])OCC[N+](C)(C)C)(COC)OC(=O)
    """)
    if core_pattern is None:
        return False, "Invalid core SMARTS pattern"
    
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Core structure not found"
    
    # Verify sn-1 ether is an alkyl chain (no adjacent carbonyl)
    ether_oxygen = None
    for match in mol.GetSubstructMatches(Chem.MolFromSmarts("[C@H](COP(=O)([O-])OCC[N+](C)(C)C)(CO[OX2])OC(=O)")):
        ether_oxygen = match[2]  # Assuming the oxygen in COC is at index 2
        break
    if ether_oxygen is None:
        return False, "Ether oxygen not found"
    
    # Check that ether oxygen is connected to a carbon chain without carbonyl
    oxygen_atom = mol.GetAtomWithIdx(ether_oxygen)
    for neighbor in oxygen_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and not any(a.GetAtomicNum() == 8 and a.GetTotalNumHs() == 0 for a in neighbor.GetNeighbors()):
            continue  # Alkyl chain
        else:
            return False, "sn-1 oxygen not connected to alkyl chain"
    
    # Verify sn-2 ester is an acyl chain (carbonyl adjacent)
    ester_oxygen = None
    for match in mol.GetSubstructMatches(Chem.MolFromSmarts("[C@H](COP(=O)([O-])OCC[N+](C)(C)C)(COC)OC(=O)")):
        ester_oxygen = match[-1]  # Assuming the oxygen in OC(=O) is last
        break
    if ester_oxygen is None:
        return False, "Ester oxygen not found"
    
    ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen)
    carbonyl_carbon = None
    for neighbor in ester_oxygen_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and any(bond.GetBondType() == Chem.BondType.DOUBLE for bond in neighbor.GetBonds()):
            carbonyl_carbon = neighbor
            break
    if carbonyl_carbon is None:
        return False, "Ester group missing carbonyl"
    
    # Check molecular weight (typical examples are >600)
    mol_wt = AllChem.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low"
    
    return True, "Has 1-alkyl-2-acyl-sn-glycero-3-phosphocholine structure"