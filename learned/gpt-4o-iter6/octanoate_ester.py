"""
Classifies: CHEBI:87657 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester features the presence of the octanoyl group `(C(=O)OCCCCCCCC)`.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for octanoyl ester linkage (-C(=O)OCCCCCCCC)
    octanoyl_ester_pattern = Chem.MolFromSmarts("C(=O)OCCCCCCCC")
    matches = mol.GetSubstructMatches(octanoyl_ester_pattern)
    
    if len(matches) == 0:
        return False, "No octanoyl ester linkage found"
    
    # Further verification of the carbon chain length
    for match in matches:
        carbon_chain_atoms = match[3:]  # The last eight atoms should be the carbon chain
        carbon_chain_length = 0
        for atom_idx in carbon_chain_atoms:
            if mol.GetAtomWithIdx(atom_idx).GetAtomicNum() == 6:  # Check for carbon atoms
                carbon_chain_length += 1
        
        if carbon_chain_length == 8:
            return True, "Contains valid octanoyl ester linkage"
    
    return False, "Octanoyl group found, but carbon chain length is incorrect"