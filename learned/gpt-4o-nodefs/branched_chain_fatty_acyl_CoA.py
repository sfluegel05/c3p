"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Attempt to parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A (CoA) substructure pattern
    coa_pattern = Chem.MolFromSmarts('C(=O)SCCNC(=O)CCNC(=O)[C@H](O)[C@H](COP(O)(=O)OP(O)(=O)O)C')

    # Check if the molecule contains Coenzyme A (CoA) moiety
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Lacks CoA substructure"

    # Possible branched chains - detect methyl branching off main chain
    branch_pattern = Chem.MolFromSmarts('CC(C)C')
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No branched chain detected"

    # If criteria are met, return True
    return True, "Molecule contains CoA moiety and branched-chain structure"