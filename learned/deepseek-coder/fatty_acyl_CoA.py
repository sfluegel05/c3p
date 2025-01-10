"""
Classifies: CHEBI:37554 fatty acyl-CoA
"""
"""
Classifies: CHEBI:15525 fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a fatty acyl-CoA based on its SMILES string.
    A fatty acyl-CoA consists of a fatty acid chain attached to coenzyme A via a thioester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a thioester bond (C(=O)S)
    thioester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[SX2]")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond found"

    # Check for the presence of the CoA moiety (adenine, ribose, and phosphate groups)
    # This pattern matches the core CoA structure including the adenine, ribose, and phosphate groups
    coa_pattern = Chem.MolFromSmarts("N1C=NC2=C1N=CN2[C@H]1O[C@H](COP(=O)(O)OP(=O)(O)O[C@H]2O[C@H](CO)[C@@H](O)[C@H]2O)[C@@H](O)[C@H]1O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Check for fatty acid chain characteristics
    # Look for at least 4 consecutive carbons in a chain
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4]")
    if not mol.HasSubstructMatch(fatty_acid_pattern):
        return False, "No fatty acid chain found"

    # Check molecular weight - fatty acyl-CoAs are typically >700 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 700:
        return False, "Molecular weight too low for fatty acyl-CoA"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbons for fatty acyl-CoA"
    if o_count < 10:
        return False, "Too few oxygens for fatty acyl-CoA"

    return True, "Contains thioester bond, CoA moiety, and fatty acid chain"