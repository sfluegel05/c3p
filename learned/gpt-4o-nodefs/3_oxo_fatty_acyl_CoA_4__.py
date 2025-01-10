"""
Classifies: CHEBI:57347 3-oxo-fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a 3-oxo-fatty acyl-CoA(4-) based on its SMILES string.
    This class includes molecules with a 3-oxo group, a fatty acyl chain, and a CoA moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3-oxo-fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the 3-oxo group pattern (C=O in a carbon chain context)
    oxo_pattern = Chem.MolFromSmarts("C=O")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not found"
    
    # Check for a fatty acyl chain (long hydrocarbon chain)
    chain_pattern = Chem.MolFromSmarts("C(C)(C)CCC")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 6:  # Typically, fatty acyl chains have at least 6 matching substructures
        return False, "Fatty acyl chain too short"
    
    # Check for CoA signature including critical components
    coa_pattern = Chem.MolFromSmarts("NC1=NC=CN=C1[C@@H]2O[C@@H]([C@H](O)[C@H]2O)COP(=O)(O)OC[C@H]3O[C@H]([C@@H](O)[C@H]3O)OP(=O)(O)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA structure not found"
    
    return True, "Contains 3-oxo group, fatty acyl chain, and CoA moiety"