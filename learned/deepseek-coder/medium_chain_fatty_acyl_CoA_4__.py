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

    # Check for CoA backbone pattern
    coa_pattern = Chem.MolFromSmarts("[O-]P(=O)([O-])OP(=O)([O-])OCC1OC(C(O)C1OP(=O)([O-])[O-])n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA backbone found"

    # Check for thioester linkage (fatty acid to CoA)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=O)[SX2]CCNC(=O)CCNC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Count carbons in fatty acid chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Medium-chain fatty acids typically have 6-12 carbons
    if not (6 <= carbon_count <= 12):
        return False, f"Carbon count {carbon_count} not in medium-chain range (6-12)"

    # Check for 4 negative charges (deprotonated phosphate groups)
    negative_charges = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1)
    if negative_charges != 4:
        return False, f"Found {negative_charges} negative charges, need exactly 4"

    # Check molecular weight - should be in expected range for medium-chain CoA
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (700 <= mol_wt <= 900):
        return False, f"Molecular weight {mol_wt:.1f} outside expected range (700-900)"

    return True, "Contains medium-chain fatty acid (6-12 carbons) linked to CoA with 4 negative charges"