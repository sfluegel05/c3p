"""
Classifies: CHEBI:61905 short-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:35601 short-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_short_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acyl-CoA based on its SMILES string.

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
    
    # Look for CoA moiety pattern
    coa_pattern = Chem.MolFromSmarts("[N&R]1(C2=NC=NC2=N)C(C(C(O1)COP(=O)(O)OP(=O)(O)OCC3OC(N4C=NC5=C4N=CN=C5N)C(O)C3O)O)O")
    coa_matches = mol.GetSubstructMatches(coa_pattern)
    if not coa_matches:
        return False, "No CoA moiety found"
    
    # Look for fatty acid moiety pattern (carbon chain with functional group)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])C[CX4][CX4][CX4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if not fatty_acid_matches:
        return False, "No fatty acid moiety found"
    
    # Check if fatty acid moiety is connected to CoA via thioester bond
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) != 1:
        return False, f"Found {len(thioester_matches)} thioester bonds, expected 1"
    
    # Count carbon atoms in the fatty acid chain
    chain_length = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C" and atom.GetIsInRingSize() == 0)
    if chain_length < 3 or chain_length > 8:
        return False, "Fatty acid chain length not in the range of 3-8 carbons"
    
    # Check for common functional groups on the fatty acid chain
    functional_groups = {"[OH]", "[NH2]", "[CH3][CX3](=[OX1])"}
    has_functional_group = any(mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)) for smarts in functional_groups)
    if not has_functional_group:
        return False, "No common functional group found on the fatty acid chain"
    
    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700 or mol_wt > 1200:
        return False, "Molecular weight outside the expected range for short-chain fatty acyl-CoA"
    
    return True, "Contains a CoA moiety connected to a short-chain fatty acid via a thioester bond"