"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
"""
Classifies: CHEBI:XXXXX medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

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

    # More flexible CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("[O-]P(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"

    # More flexible thioester pattern
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]CCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Identify fatty acid chain
    # Match from thioester carbon to first non-carbonyl carbon
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]CCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    fatty_acid_match = mol.GetSubstructMatch(fatty_acid_pattern)
    if not fatty_acid_match:
        return False, "Could not identify fatty acid chain"

    # Count carbons in fatty acid chain
    # Start from thioester carbon and count until CoA backbone
    start_atom = mol.GetAtomWithIdx(fatty_acid_match[0])
    carbon_count = 0
    current_atom = start_atom
    visited = set()
    
    while True:
        if current_atom.GetAtomicNum() == 6:
            carbon_count += 1
        visited.add(current_atom.GetIdx())
        
        # Get neighbors not in CoA backbone
        neighbors = [n for n in current_atom.GetNeighbors() 
                    if n.GetIdx() not in fatty_acid_match and n.GetIdx() not in visited]
        
        if not neighbors:
            break
        current_atom = neighbors[0]

    # Medium-chain fatty acids typically have 6-12 carbons
    if not (6 <= carbon_count <= 12):
        return False, f"Fatty acid chain has {carbon_count} carbons, not in medium-chain range (6-12)"

    # Check for 4 negative charges (deprotonated phosphate groups)
    negative_charges = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1)
    if negative_charges != 4:
        return False, f"Found {negative_charges} negative charges, need exactly 4"

    # More flexible molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (600 <= mol_wt <= 1000):
        return False, f"Molecular weight {mol_wt:.1f} outside expected range (600-1000)"

    return True, f"Contains medium-chain fatty acid ({carbon_count} carbons) linked to CoA with 4 negative charges"