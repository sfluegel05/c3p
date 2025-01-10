"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:75840 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is a fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"

    # Check if any ester group is attached to a 10-carbon chain
    for match in ester_matches:
        ester_oxygen = match[0]
        ester_carbon = match[1]
        
        # Get the atom connected to the ester carbon (should be part of the 10-carbon chain)
        ester_carbon_neighbors = mol.GetAtomWithIdx(ester_carbon).GetNeighbors()
        for neighbor in ester_carbon_neighbors:
            if neighbor.GetIdx() != ester_oxygen:
                # Traverse the carbon chain starting from this neighbor
                chain_length = 0
                current_atom = neighbor
                visited = set()
                while current_atom and current_atom.GetAtomicNum() == 6:  # Carbon atom
                    if current_atom.GetIdx() in visited:
                        break  # Avoid infinite loops
                    visited.add(current_atom.GetIdx())
                    chain_length += 1
                    if chain_length >= 10:
                        break
                    # Move to the next carbon in the chain
                    next_atom = None
                    for next_neighbor in current_atom.GetNeighbors():
                        if next_neighbor.GetAtomicNum() == 6 and next_neighbor.GetIdx() not in visited:
                            next_atom = next_neighbor
                            break
                    current_atom = next_atom
                
                if chain_length >= 10:
                    return True, "Contains a decanoate ester group (10-carbon chain with ester linkage)"

    return False, "No decanoate ester group found (10-carbon chain with ester linkage)"