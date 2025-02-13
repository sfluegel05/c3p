"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: CHEBI:38034 acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester that results from the formal condensation of 
    the thiol group of coenzyme A with the carboxy group of any carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for CoA backbone using multiple patterns
    coa_patterns = [
        Chem.MolFromSmarts("[N]1C=NC2=C1N=CN=C2N[C@@H]3[C@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(OCC[N]4C=NC5=C4N=CN=C5N)O)O)O"),
        Chem.MolFromSmarts("[N]1C=NC2=C1N=CN=C2N[C@@H]3[C@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)OCC[N]4C=NC5=C4N=CN=C5N)O)O"),
        Chem.MolFromSmarts("[N]1C=NC2=C1N=CN=C2N[C@@H]3[C@H]([C@@H]([C@H](O3)COP(=O)(O)OP(=O)(O)CCC[N]4C=NC5=C4N=CN=C5N)O)O")
    ]
    
    coa_match = any(mol.HasSubstructMatch(pattern) for pattern in coa_patterns)
    if not coa_match:
        return False, "No CoA backbone found"
    
    # Look for thioester (-C(=O)S-) group
    thioester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester group found"
    
    # Look for carboxylic acid (-C(=O)O-) group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxylic acid group found"
    
    return True, "Contains CoA backbone and thioester linking a carboxylic acid"