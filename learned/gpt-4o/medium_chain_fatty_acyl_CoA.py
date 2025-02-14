"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA is a fatty acyl-CoA with a fatty acid chain of 6-12 C atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Complete pattern for Coenzyme A core structure including adenine, ribose, and thiol linkage
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)[C@H](COP(=O)(O)OP(=O)(O)OC[C@@H]1O[C@@H](O)[C@@H](O)C1O)n1cnc2c(ncnc2n1)C3CC(O3)O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Missing complete Coenzyme A structure"

    # Check for thioester linkage (should be part of the fatty acid-CoA linkage): -C(=O)-S- with linkage to a medium chain
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found to CoA or incorrect fatty chain length"

    # Extract carbon chain length in fatty acyl chain
    fatty_acid_chain = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    chain_length = len(fatty_acid_chain)
    if not (6 <= chain_length <= 12):
        return False, f"Fatty acid chain length not within 6-12 carbons, found {chain_length} carbons"

    return True, "Contains medium-chain fatty acyl-CoA structure with correct components and chain length"