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

    # Define more refined CoA substructure pattern
    # Including adenine part, ribose, pyrophosphate, and the expected acyl chain linker with 'S'
    coa_pattern = Chem.MolFromSmarts('OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(=O)O)O[C@H]2n3cnc2c(N)ncn3C'
                                     'NC(=O)CCNC(=O)[C@H](O)C(C)(C)[O-]')  
    
    # Check for the presence of CoA substructure
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Lacks CoA substructure"

    # Refined branching patterns to capture various methyl groupings
    branch_patterns = [
        Chem.MolFromSmarts('CC(C)C'),       # Simple methyl branching
        Chem.MolFromSmarts('C(C)(C)C'),     # Double methyl branching
        Chem.MolFromSmarts('[CH2]C(C)C=O'), # Keto group on branched channel
        Chem.MolFromSmarts('C=C(C)C'),      # Brancing with double bond
        Chem.MolFromSmarts('CCC(C)C')       # Longer chain to capture larger structures
    ]
    
    # Check for any branched patterns
    branched_match = any(mol.HasSubstructMatch(p) for p in branch_patterns)
    if not branched_match:
        return False, "No branched chain detected"

    return True, "Molecule contains CoA moiety and branched-chain structure"