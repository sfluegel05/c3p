"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.
    A medium-chain fatty acyl-CoA is a CoA derivative where a medium-chain fatty acid (6-12 carbons)
    is attached to CoA via a thioester bond.

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

    # Look for the thioester bond (S-C(=O)-)
    thioester_pattern = Chem.MolFromSmarts("[SX2][CX3](=O)")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester bond found"

    # Look for the CoA moiety (simplified pattern)
    coa_pattern = Chem.MolFromSmarts("[SX2][CX3](=O)[CX4][NX3][CX4][CX4](=O)[NX3][CX4][CX4](=O)[OX2]")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Identify the fatty acid chain attached to the thioester bond
    # The chain should have 6-12 carbons
    chain_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) < 1:
        return False, "No fatty acid chain found"

    # Count the number of carbons in the chain
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 6 or c_count > 12:
        return False, f"Chain length {c_count} is not within the medium-chain range (6-12 carbons)"

    return True, "Contains a medium-chain fatty acid attached to CoA via a thioester bond"