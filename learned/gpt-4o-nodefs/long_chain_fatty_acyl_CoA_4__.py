"""
Classifies: CHEBI:83139 long-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_long_chain_fatty_acyl_CoA_4_(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acyl-CoA(4-) based on its SMILES string.

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

    # Flexible pattern for a long hydrophobic tail (longer chain with possible unsaturations or branching)
    carbon_chain_pattern = Chem.MolFromSmarts("C[C@H](CCCCCCCCCCCC=*)*")  # loosely matches long carbon chains
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No appropriate long carbon chain typical of fatty acid found"

    # Check for entire CoA moiety pattern, it should include thioester linkage and adenosine
    coa_full_pattern = Chem.MolFromSmarts("SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(OP(=O)(OC[C@H]1O[C@@H](N2C3=NC(N)=NC(=N3)N)[C@H]([C@@H]1O)OP(O)(O)=O)O)O")
    if not mol.HasSubstructMatch(coa_full_pattern):
        return False, "No complete CoA moiety detected"

    # Check for the presence of multiple phosphate groups, the expected pattern in CoA
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(OP(=O)([O-])[O-])[O-]")
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    if phosphate_count < 3:
        return False, f"Insufficient phosphate groups detected, found {phosphate_count}"

    return True, "Matches structure of long-chain fatty acyl-CoA(4-) with complete CoA moiety and long carbon chain"