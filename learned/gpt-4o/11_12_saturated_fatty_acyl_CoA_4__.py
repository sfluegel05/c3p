"""
Classifies: CHEBI:84948 11,12-saturated fatty acyl-CoA(4-)
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_11_12_saturated_fatty_acyl_CoA_4__(smiles):
    """
    Determines if a molecule is an 11,12-saturated fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a corrected and comprehensive CoA SMARTS pattern
    coa_core_smarts = "SCCNC(=O)CCNC(=O)[C@H](O)C(C)COP(=O)(O)OP(=O)(O)OCC1OC(CO)C(O)C1O"
    coa_pattern = Chem.MolFromSmarts(coa_core_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A group not found"

    # Define a general pattern for carbon acyl chains to detect the 11th to 12th saturation
    # This pattern should ensure a single bond between the 11th and 12th carbon
    general_chain_pattern = Chem.MolFromSmarts('CCCCCCCCCCCCCCCCCC')
    chain_matches = mol.GetSubstructMatches(general_chain_pattern)

    if not chain_matches:
        return False, "No valid long carbon chain found"

    for match in chain_matches:
        # Check 11th to 12th C-C bond for saturation
        if len(match) >= 12:
            bond = mol.GetBondBetweenAtoms(match[10], match[11])
            if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                return True, "Contains valid long carbon chain with 11-12 bond saturation within fatty acyl part"

    return False, "No saturated bond between 11th and 12th carbon atoms found in carbon chain"

__metadata__ = { 'chemical_class': { 'id': 'CHEBI:xxxx', 'name': '11,12-saturated fatty acyl-CoA(4-)', 'definition': 'Any fatty acyl-CoA(4-) in which the 11-12 bond of the fatty acyl group is saturated' } }