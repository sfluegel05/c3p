"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
from rdkit import Chem

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphinganine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define sphinganine backbone pattern
    # Two adjacent hydroxyl groups and one amine group, ensuring correct chirality
    sphinganine_backbone = Chem.MolFromSmarts("[C@@H](O)[C@H](CO)N")
    if not mol.HasSubstructMatch(sphinganine_backbone):
        return False, "No sphinganine backbone found"
    
    # Define N-acyl linkage pattern with aliphatic chain constraint
    # Aliphatic chain typically contains 12-34 carbons
    acyl_linkage_pattern = Chem.MolFromSmarts("C(=O)N[C@@H](CO)C")
    acyl_matches = mol.GetSubstructMatches(acyl_linkage_pattern)
    if not acyl_matches:
        return False, "No N-acyl linkage found"
    
    # Check the length of the carbon chain
    # Assume chain should be between 12 to 34 carbons as typical N-acylsphinganines
    chain_length = 0
    for match in acyl_matches:
        atom_indices = [atom.GetIdx() for atom in mol.GetAtoms()]
        chain = [mol.GetAtomWithIdx(idx) for idx in atom_indices if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6]
        chain_length = len(chain)
        if 12 <= chain_length <= 34:
            break
    else:
        return False, f"Aliphatic chain is not in the expected length range; found {chain_length} carbons"

    # Allow optional sugars or additional groups
    possible_headgroup = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](O)C[C@H]1O")  # Simplified sugar example
    
    if mol.HasSubstructMatch(possible_headgroup):
        return True, "Contains sphinganine backbone with N-acyl linkage and possible sugar headgroup"

    return True, "Contains sphinganine backbone with N-acyl linkage"