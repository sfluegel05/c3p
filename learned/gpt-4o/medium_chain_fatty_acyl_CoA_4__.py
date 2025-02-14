"""
Classifies: CHEBI:90546 medium-chain fatty acyl-CoA(4-)
"""
from rdkit import Chem

def is_medium_chain_fatty_acyl_CoA_4__(smiles: str):
    """
    Determines if a molecule is a medium-chain fatty acyl-CoA(4-) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a medium-chain fatty acyl-CoA(4-), False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for CoA moiety
    # Updated pattern to better capture CoA structure
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # Check for the presence of thioester linkage
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Identify and examine the carbon chain length using a pattern
    chain_pattern = Chem.MolFromSmarts("C(=O)SCCCCCCCC")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    chain_lengths = [len(match) for match in chain_matches]

    if not any(chain_length >= 8 and chain_length <= 12 for chain_length in chain_lengths):
        return False, f"Chain lengths {chain_lengths} do not include a medium-chain length (8-12 carbons)"
    
    # Check that there are enough deprotonated oxygens (4- charge)
    deprotonated_oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if deprotonated_oxygen_count < 4:
        return False, f"Insufficient number of deprotonated oxygens: {deprotonated_oxygen_count}"

    return True, "Molecule is a medium-chain fatty acyl-CoA(4-)"