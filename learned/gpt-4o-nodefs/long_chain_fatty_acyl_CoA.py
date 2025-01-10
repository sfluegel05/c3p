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

    # Create a pattern for Coenzyme A-like structure, which may include nucleotide features
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)C")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing Coenzyme A-like structure"

    # Check for the presence of a long-chain fatty acid-like alkyl chain
    # This is a longer chain pattern with CH2 groups
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCCCCC")
    if not mol.HasSubstructMatch(long_chain_pattern):
        return False, "No long hydrocarbon chain found"

    # Check for the presence of thiol ester linkage
    thiol_ester_pattern = Chem.MolFromSmarts("C(=O)S")
    if not mol.HasSubstructMatch(thiol_ester_pattern):
        return False, "Missing thiol ester linkage"

    # Count total number of carbons, which should be typical of long-chain fatty acyl groups
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14:
        return False, "Too few carbons for long-chain fatty acyl group"

    return True, "Matches long-chain fatty acyl-CoA features"