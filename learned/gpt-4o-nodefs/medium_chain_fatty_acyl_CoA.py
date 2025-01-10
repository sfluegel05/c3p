"""
Classifies: CHEBI:61907 medium-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_medium_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA based on its SMILES string.

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

    # Improved Coenzyme A moiety pattern
    # A wide-capturing SMARTS for the CoA moiety
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)([O-])=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"

    # Thioester bond pattern
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Get Thioester Group and Attached Acyl Chain
    # Consider all carbon chains attached to the C=O-S group
    thioester_matches = mol.GetSubstructMatches(thioester_pattern)
    
    # Identifying Acyl Chains
    acyl_chain_lengths = []
    for match in thioester_matches:
        atom_idx = match[0]  # Starting from the carbon in C(=O)
        acyl_chain = Chem.FragmentOnBonds(mol, [atom_idx])
        acyl_atoms = [atom for atom in acyl_chain.GetAtoms() if atom.GetAtomicNum() == 6]
        acyl_length = len(acyl_atoms)
        acyl_chain_lengths.append(acyl_length)

    # Check if there is any medium chain (6-12 carbons) acyl component
    if any(6 <= length <= 12 for length in acyl_chain_lengths):
        return True, "Contains a CoA moiety with medium-chain fatty acyl group"

    return False, "Fatty acyl chain not of medium length"