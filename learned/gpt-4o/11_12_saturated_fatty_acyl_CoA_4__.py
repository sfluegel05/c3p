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

    # Define simplified CoA pattern based on known core scaffold
    coa_core_smarts = "SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(CO)C(O)C1O"
    coa_pattern = Chem.MolFromSmarts(coa_core_smarts)
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "Coenzyme A group not found"

    # Define general pattern for long carbon chains
    # Searching concretely for saturated bonds between 11th and 12th carbon
    acyl_chain_pattern = Chem.MolFromSmarts('CCCC[CH2][CH2]CCCCCCCCCCCCCC(=O)')
    chain_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    
    # Check saturation state specifically between 11th and 12th carbon
    if not chain_matches:
        return False, "No valid fatty acyl chain with appropriate saturation found"
    for match in chain_matches:
        # Extract indices of interest (11th and 12th)
        if len(match) < 12:
            continue
        # Verify saturation specifically between match[10] and match[11]
        bond = mol.GetBondBetweenAtoms(match[10], match[11])
        if bond and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            return True, "Contains CoA group with a saturated bond from 11-12 in fatty acyl chain"

    return False, "11-12 carbon bond is not saturated or chain length inadequate"

# The program remakes patterns to improve likelihood of correct recognition and detail tracing is done systematically