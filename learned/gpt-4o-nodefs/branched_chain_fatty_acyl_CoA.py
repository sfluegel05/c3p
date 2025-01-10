"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Enhanced CoA substructure pattern
    coa_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)CCNC(=O)[C@H](O)[C@@H](COP(O)(=O)OP(O)(=O)O)C')

    # Check for the presence of CoA substructure
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Lacks CoA substructure"

    # Enhanced branched chain pattern (multiple possible branching patterns)
    branch_patterns = [
        Chem.MolFromSmarts('CC(C)C'),       # Simple methyl branching
        Chem.MolFromSmarts('C(C)C(C)C'),    # Multi branched
        Chem.MolFromSmarts('[C@H](C)(C)'),  # Chiral branching
        Chem.MolFromSmarts('C=C(C)C')       # Brancing with double bond
    ]
    
    # Check for any branched patterns
    branched_match = any(mol.HasSubstructMatch(p) for p in branch_patterns)
    if not branched_match:
        return False, "No branched chain detected"

    # Confirm required functional groups or structural features
    # Check for possible hydroxylation or other characteristic groups
    hydroxyl_pattern = Chem.MolFromSmarts('[CX4](O)')  # Hydroxylation
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Potential missing hydroxyl group for characteristic feature"

    return True, "Molecule contains CoA moiety and branched-chain structure"