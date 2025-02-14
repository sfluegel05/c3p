"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
"""
Classifies: CHEBI:60444 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    A branched-chain fatty acyl-CoA results from the formal condensation of the thiol group of
    coenzyme A with the carboxy group of any branched-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone
    coa_pattern = Chem.MolFromSmarts("C(C)(COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12)C(=O)NCCCC(=O)NCCSC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No coenzyme A backbone found"
    
    # Look for branched fatty acid chain
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3][CX4,CX3]([CX4,CX3])([CX4,CX3])[CX3,CX2]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "No branched fatty acid chain found"
    
    # Check for ester linkage between CoA and fatty acid
    ester_pattern = Chem.MolFromSmarts("CCC(=O)OC")
    if not mol.HasSubstructMatch(ester_pattern):
        # Try a more relaxed ester pattern
        ester_pattern = Chem.MolFromSmarts("C(=O)OC")
        if not mol.HasSubstructMatch(ester_pattern):
            return False, "No ester linkage found"
    
    # Check for specific branched fatty acid groups
    branched_groups = ["CC(C)", "CC(C)C", "CCC(C)", "CCCC(C)", "CCCCC(C)", "CCCCCC(C)"]
    branched_group_found = False
    for group in branched_groups:
        group_pattern = Chem.MolFromSmarts(group)
        if mol.HasSubstructMatch(group_pattern):
            branched_group_found = True
            break
    if not branched_group_found:
        return False, "No specific branched fatty acid group found"
    
    # Exclude molecules with hydroxyl groups on the fatty acid chain
    hydroxyl_pattern = Chem.MolFromSmarts("CC(O)")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Hydroxyl group present on the fatty acid chain"
    
    # Check molecular weight - typically >800 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 800:
        return False, "Molecular weight too low for branched-chain fatty acyl-CoA"
    
    return True, "Contains coenzyme A backbone with a branched-chain fatty acid attached via an ester bond"