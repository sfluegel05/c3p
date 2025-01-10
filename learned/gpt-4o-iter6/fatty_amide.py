"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is defined by an amide group derived from a long-chain fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of amide bond (C(=O)N)
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide bond found"

    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Check for long carbon chain after the carbonyl group
    for match in amide_matches:
        carbonyl_c_idx = match[0]
        nitrogen_idx = match[1]
        
        # Start from the carbonyl carbon and assess chain length connected to nitrogen
        visited = set()
        carbon_chain_length = 0

        # Check carbon atoms connected to carbonyl carbon, excluding the direct amide nitrogen
        for atom in mol.GetAtomWithIdx(carbonyl_c_idx).GetNeighbors():
            if atom.GetIdx() != nitrogen_idx:
                queue = [(atom, 1)]  # (atom, depth)
        
                while queue:
                    current_atom, depth = queue.pop(0)
                    if current_atom.GetIdx() in visited or depth > 2:
                        continue
                    
                    visited.add(current_atom.GetIdx())

                    if current_atom.GetAtomicNum() == 6:  # Check if it's a carbon atom
                        carbon_chain_length += 1

                    # Only add carbon atoms (ignores branching through heteroatoms)
                    queue.extend([(neighbor, depth + 1) for neighbor in current_atom.GetNeighbors() if neighbor.GetIdx() not in visited and neighbor.GetAtomicNum() == 6])

        # Verify if enough carbon atoms are found
        if carbon_chain_length >= 8:
            return True, "Contains amide bond derived from a long chain fatty acid"

    return False, "Insufficient carbon chain length beyond amide bond"