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
    
    # Updated CoA pattern to capture key structural motifs more robustly
    coa_pattern = Chem.MolFromSmarts("NC(=O)CCNC(=O)[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H](n2cnc3c(N)ncnc23)[C@H](O)[C@@H]1OP([O-])([O-])=O")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No CoA moiety found"
    
    # Check for the presence of thioester linkage
    # Expanded pattern to ensure coverage
    thioester_pattern = Chem.MolFromSmarts("C(=O)SCCC")
    if not mol.HasSubstructMatch(thioester_pattern):
        return False, "No thioester linkage found"

    # Check for carbon chain length: look for C(=O)S followed by the carbon chain
    chains = mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)SCC"))
    chain_lengths = []
    
    for match in chains:
        chain_length = 0
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetSymbol() == 'C':
                chain_length += 1
        # Account the carbon atoms in the backbone after the thioester
        chain_lengths.append(chain_length)

    if not any(chain_length >= 8 and chain_length <= 12 for chain_length in chain_lengths):
        return False, f"Chain lengths {chain_lengths} do not include a medium-chain length (8-12 carbons)"
    
    # Check that there are enough deprotonated oxygens (4- charge)
    deprotonated_oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1)
    if deprotonated_oxygen_count < 4:
        return False, f"Insufficient number of deprotonated oxygens: {deprotonated_oxygen_count}"

    return True, "Molecule is a medium-chain fatty acyl-CoA(4-)"