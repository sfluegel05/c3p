"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.
    Long-chain fatty acyl-CoA(4-) are characterized by a long fatty acid chain attached to a 
    Coenzyme A structure with a net charge of 4- due to deprotonated phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define an improved pattern for Coenzyme A
    coa_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]([C@@H]1OP([O-])=O)n1cnc2c(N)ncnc12")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A structure found or incorrect stereochemistry"

    # Define a more flexible pattern for long aliphatic chains, considering unsaturations
    fatty_acyl_pattern = Chem.MolFromSmarts("C(=O)SC[CH2;R0][CX4,CX3;H1,H2]~[CX4,CX3]{10,}")
    if not mol.HasSubstructMatch(fatty_acyl_pattern):
        return False, "No long-chain fatty acyl group detected"

    # Ensure the presence of deprotonated phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("[P](=O)([O-])([O-])[O]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 3:
        return False, "Insufficient number of deprotonated phosphate groups"

    return True, "Classified as a long-chain fatty acyl-CoA(4-)"