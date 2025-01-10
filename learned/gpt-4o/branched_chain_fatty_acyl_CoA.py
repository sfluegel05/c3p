"""
Classifies: CHEBI:61912 branched-chain fatty acyl-CoA
"""
from rdkit import Chem

def is_branched_chain_fatty_acyl_CoA(smiles: str):
    """
    Classifies if a molecule is a branched-chain fatty acyl-CoA based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Coenzyme A structure pattern including pantetheine and nucleotide parts
    coa_pattern = Chem.MolFromSmarts("C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(=O)(O)O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A structure not found"
    
    # Thioester bond pattern ensuring connection to CoA
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCNC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "Thioester linkage associated with CoA not found"

    # Improved branching detection: Look for isopropyl/iso-butyl groups on aliphatic chains
    branch_pattern = Chem.MolFromSmarts("C(C)(C)[CH2,CH,CH3]")
    if not mol.HasSubstructMatch(branch_pattern):
        return False, "No branched alkyl chain found"

    return True, "Contains branched-chain fatty acyl-CoA structure"