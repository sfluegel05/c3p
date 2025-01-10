"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine contains a six-membered ring with two nitrogens and 
    two carbonyl groups at the 2 and 5 positions with various substitutions possible.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create SMARTS pattern for 2,5-diketopiperazine
    # The pattern ensures a six-membered ring with two adjacent carbonyl groups connected to two nitrogens
    diketopiperazine_pattern = Chem.MolFromSmarts('C1C(N)C(=O)NC(=O)C1')  # Explicit 6-membered ring

    # Check ring substructure match
    if not mol.HasSubstructMatch(diketopiperazine_pattern):
        return False, "Does not contain 2,5-diketopiperazine skeleton"

    # Additional checks: Ensure matching ring size and connectivity
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) == 6:
            atom_types = [mol.GetAtomWithIdx(idx).GetSymbol() for idx in ring]
            # Should include two nitrogens and at least two oxygens for diketopiperazine
            if atom_types.count('N') == 2 and atom_types.count('O') >= 2:
                return True, "Contains valid 2,5-diketopiperazine skeleton"
    
    return False, "Ring does not match 2,5-diketopiperazine configuration"