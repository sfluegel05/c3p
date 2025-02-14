"""
Classifies: CHEBI:139575 monounsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monounsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acyl-CoA based on its SMILES string.
    A monounsaturated fatty acyl-CoA has one carbon-carbon double bond in the fatty acid chain, and a CoA moiety linked via a thioester.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the CoA moiety using a simplified SMARTS focusing on the core.
    coa_pattern = Chem.MolFromSmarts("[P](=O)([O])O[C@H]1[C@@H]([C@@H]([C@H](O1)C2=CN=C3N=CN=C(N)C3=C2)O)COP(=O)([O])O") #Simplified CoA pattern
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA moiety not found"

    # 2. Check for the thioester bond (S-C=O) - this should be connected to the acyl chain
    thioester_pattern = Chem.MolFromSmarts("SC(=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester bond not found"

    # 3. Check for exactly one carbon-carbon double bond
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_matches = mol.GetSubstructMatches(double_bond_pattern)
    if len(double_bond_matches) != 1:
        return False, f"Molecule has {len(double_bond_matches)} carbon-carbon double bonds, expected 1"
    
    # 4. Check for a long carbon chain
    fatty_acid_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_chain_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing a fatty acid chain"

    # 5. Check for number of rotatable bonds - should be relatively long
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 8:
        return False, "Fatty acid chain too short"


    return True, "Molecule is a monounsaturated fatty acyl-CoA"