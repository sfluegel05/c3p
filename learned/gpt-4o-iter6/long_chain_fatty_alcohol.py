"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol has a carbon chain length ranging from C13 to C22 with
    at least one hydroxyl (-OH) group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for hydroxyl group (-OH) presence
    oh_group = Chem.MolFromSmarts("[CX4][OX2H]")
    if not mol.HasSubstructMatch(oh_group):
        return False, "No hydroxyl group found"
    
    # Identify the longest chain of carbon atoms
    carbon_chains = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        bonds = [b for b in atom.GetBonds() if b.GetOtherAtom(atom).GetAtomicNum() == 6]
        chain_length = len(bonds) + 1
        carbon_chains.append(chain_length)
    
    # Check if any chain meets the C13 to C22 requirement
    max_chain_length = max(carbon_chains) if carbon_chains else 0
    if 13 <= max_chain_length <= 22:
        return True, f"Contains a carbon chain of length {max_chain_length} and a hydroxyl group"
    else:
        return False, f"Carbon chain length of {max_chain_length}, required between 13 and 22"