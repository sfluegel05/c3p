"""
Classifies: CHEBI:35915 sterol ester
"""
from rdkit import Chem

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined Sterol core identification: a core structure invoking common sterol patterns.
    # Sterols have a tetracyclic skeleton with hydrocarbon chains and a hydroxyl group typically on C3.
    sterol_smarts = "C1CCC2C3C(CC(C4C3(CCC2C1)C)C)CC4" # General sterol scaffold
    sterol_core = Chem.MolFromSmarts(sterol_smarts)

    if not mol.HasSubstructMatch(sterol_core):
        return False, "No sterol core structure found"

    # Ester group pattern (-C(=O)O-)
    ester_pattern_smarts = "C(=O)O" # Generic ester linkage
    ester_pattern = Chem.MolFromSmarts(ester_pattern_smarts)
    
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Count carbons in ester-linked aliphatic chain to ensure fatty acid presence
    for match in mol.GetSubstructMatches(ester_pattern):
        ester_atom_idx = match[1]  # The oxygen in the ester linkage
        
        # Follow the carbon chain from here:
        carbon_chain_length = 0
        atom = mol.GetAtomWithIdx(ester_atom_idx)
        
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon neighbor
                carbon_count = recursive_carbon_count(neighbor, {atom.GetIdx()})
                if carbon_count >= 10:  # Typical minimum length for fatty acids
                    return True, "Contains sterol core structure esterified with a sufficiently long aliphatic chain"

    return False, "No sufficiently long aliphatic chain found"

def recursive_carbon_count(atom, visited):
    """
    Helper function to recursively count the number of carbons in a chain.
    """
    count = 0
    
    # Traverse through neighbors
    for neighbor in atom.GetNeighbors():
        if neighbor.GetIdx() not in visited:
            visited.add(neighbor.GetIdx())
            if neighbor.GetAtomicNum() == 6:  # Carbon
                count += 1 + recursive_carbon_count(neighbor, visited)
    
    return count

# Example metadata
__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:None',
        'name': 'sterol ester',
        'definition': 'Chemical entities consisting of a sterol esterified with a fatty acid'
    },
    'config': {
        'llm_model_name': 'rdkit_based_model',
        'f1_threshold': 0.8,
        'max_attempts': 5
    },
    'message': None,
    'attempt': 3,
    'success': True
}