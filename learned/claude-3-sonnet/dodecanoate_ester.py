"""
Classifies: CHEBI:87659 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from collections import deque

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester contains a lauric acid (C12:0) component attached via an ester bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    def count_carbon_chain(start_atom, visited=None):
        """Helper function to count continuous carbon chain length from a starting atom"""
        if visited is None:
            visited = set()
            
        queue = deque([(start_atom, 1)])  # (atom, chain_length)
        max_length = 1
        visited.add(start_atom.GetIdx())
        
        while queue:
            current_atom, length = queue.popleft()
            max_length = max(max_length, length)
            
            for neighbor in current_atom.GetNeighbors():
                if (neighbor.GetAtomicNum() == 6 and  # is carbon
                    neighbor.GetIdx() not in visited and
                    sum(1 for bond in neighbor.GetBonds() if bond.GetBondType() != Chem.BondType.SINGLE) == 0):  # saturated
                    visited.add(neighbor.GetIdx())
                    queue.append((neighbor, length + 1))
                    
        return max_length

    for match in ester_matches:
        carbonyl_carbon = mol.GetAtomWithIdx(match[1])  # Get the C(=O) carbon
        
        # Check each carbon attached to the carbonyl carbon
        for neighbor in carbonyl_carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # is carbon
                # Count the length of the carbon chain
                chain_length = count_carbon_chain(neighbor)
                
                if chain_length >= 12:
                    # Verify it's a saturated chain by checking for any double bonds
                    atoms_to_check = {atom.GetIdx() for atom in mol.GetAtoms() 
                                   if atom.GetAtomicNum() == 6}  # all carbons
                    has_double_bond = False
                    
                    for bond in mol.GetBonds():
                        if (bond.GetBondType() == Chem.BondType.DOUBLE and
                            bond.GetBeginAtomIdx() in atoms_to_check and
                            bond.GetEndAtomIdx() in atoms_to_check):
                            has_double_bond = True
                            break
                    
                    if not has_double_bond:
                        return True, "Contains dodecanoate (laurate) ester group with saturated C12 chain"
    
    return False, "No dodecanoate ester group with appropriate chain length found"