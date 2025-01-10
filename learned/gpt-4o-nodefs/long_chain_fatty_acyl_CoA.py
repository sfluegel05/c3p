"""
Classifies: CHEBI:33184 long-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_long_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create a pattern for Coenzyme A-like structure
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A-like structure"

    # Check for the presence of a thiol ester linkage and COA pattern
    thiol_ester_coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)")
    if not mol.HasSubstructMatch(thiol_ester_coa_pattern):
        return False, "Missing proper thiol ester and CoA linkage"

    # Check for the presence of a long-chain (12-24 carbons) fatty acid-like alkyl chain
    long_chain_pattern = Chem.MolFromSmarts("C(C)(C)COP([O-])([O-])=O")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Check number of double bonds in the chain
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bond_count = len(mol.GetSubstructMatches(double_bond_pattern))
    if double_bond_count > 5:
        return False, "Too many double bonds for typical long-chain fatty acyl-CoA"

    # Count total number of carbons, which should be typical of long fatty acyl groups
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not 16 <= c_count <= 24:
        return False, "Carbon count not typical for long-chain fatty acyl group"

    return True, "Matches long-chain fatty acyl-CoA features"