"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    A long-chain fatty acyl-CoA(4-) is a molecule with a long-chain fatty acid (14-24 carbons)
    attached to a CoA moiety via a thioester bond, with deprotonated phosphate and diphosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for CoA moiety pattern (deprotonated phosphate and diphosphate groups)
    coa_pattern = Chem.MolFromSmarts("[O-]P(=O)([O-])OP(=O)([O-])OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety with deprotonated phosphate and diphosphate groups found"

    # Check for thioester bond (S-C(=O)-)
    thioester_pattern = Chem.MolFromSmarts("[SX2][CX3](=[OX1])")
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    if len(thioester_matches) == 0:
        return False, "No thioester bond found"

    # Check for long-chain fatty acid (14-24 carbons)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) == 0:
        return False, "No long-chain fatty acid found"

    # Count carbons in the fatty acid chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14 or c_count > 24:
        return False, f"Fatty acid chain length ({c_count} carbons) is not within the expected range (14-24 carbons)"

    # Check molecular weight - long-chain fatty acyl-CoA(4-) typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for long-chain fatty acyl-CoA(4-)"

    return True, "Contains a long-chain fatty acid attached to a CoA moiety via a thioester bond, with deprotonated phosphate and diphosphate groups"