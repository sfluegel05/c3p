"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broadened CoA substructure pattern to include more specific phosphate and base features
    coa_pattern = Chem.MolFromSmarts('NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC')

    # Check for the presence of CoA substructure
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Lacks CoA substructure"

    # Expanded branching patterns to capture various methyl groupings
    branch_patterns = [
        Chem.MolFromSmarts('CC(C)C'),       # Simple methyl branching
        Chem.MolFromSmarts('C(C)(C)C'),     # Double methyl branching
        Chem.MolFromSmarts('CCC(C)C'),      # Additional chiral branching
        Chem.MolFromSmarts('C=C(C)C'),      # Brancing with double bond
        Chem.MolFromSmarts('CC(O)C')        # Hydroxyl attached to branch point
    ]
    
    # Check for any branched patterns
    branched_match = any(mol.HasSubstructMatch(p) for p in branch_patterns)
    if not branched_match:
        return False, "No branched chain detected"

    # It might or might not have hydroxylation, broaden functional group checks
    return True, "Molecule contains CoA moiety and branched-chain structure"