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
    
    # Look for the 3-oxo group pattern (C=O as part of a carbon chain)
    oxo_pattern = Chem.MolFromSmarts("CC(=O)C")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "3-oxo group not found"
    
    # Check for a fatty acyl chain (long hydrocarbon chain)
    chain_pattern = Chem.MolFromSmarts("CCCC")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) <= 6: # typically fatty acyl chains have at least 6 carbons
        return False, "Fatty acyl chain too short"
    
    # Check for CoA signature (including adenine and phosphates)
    coa_pattern = Chem.MolFromSmarts("NC1=C2N=CN([C@H]3O[C@H]([C@H](O)[C@@H]3OP(=O)([O-])O)COP(=O)([O-])O)C2=NC=N1")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "CoA structure not found"

    return True, "Contains 3-oxo group, fatty acyl chain, and CoA moiety"