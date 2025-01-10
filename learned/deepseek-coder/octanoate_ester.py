"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:75548 octanoate ester
"""
from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).

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

    # Define the octanoate ester pattern: C(=O)O-R where R is any group
    octanoate_pattern = Chem.MolFromSmarts("CCCCCCCC(=O)O[!H]")
    
    # Check if the pattern is present in the molecule
    if not mol.HasSubstructMatch(octanoate_pattern):
        return False, "No octanoate ester group found"

    # Ensure the ester group is derived from octanoic acid
    # Check that the carbon chain is exactly 8 carbons long
    ester_matches = mol.GetSubstructMatches(octanoate_pattern)
    for match in ester_matches:
        ester_atom = mol.GetAtomWithIdx(match[0])  # Carbon in C=O
        if ester_atom.GetTotalDegree() != 3:  # Ensure it's part of an ester
            return False, "Not a valid ester group"
        
        # Check the chain length
        chain_length = 0
        neighbor = ester_atom.GetNeighbors()[0]  # First neighbor is the chain
        while neighbor.GetAtomicNum() == 6:  # Carbon
            chain_length += 1
            if chain_length > 7:  # 7 carbons + 1 in C=O = 8 carbons total
                break
            # Move to next carbon, ensuring it's part of a linear chain
            next_neighbors = [n for n in neighbor.GetNeighbors() if n.GetIdx() != ester_atom.GetIdx()]
            if len(next_neighbors) != 1:
                return False, "Chain is branched or not linear"
            neighbor = next_neighbors[0]
        
        if chain_length != 7:  # 7 carbons + 1 in C=O = 8 carbons total
            return False, "Chain length not consistent with octanoic acid"

    return True, "Contains octanoate ester group (CCCCCCCC(=O)OR)"