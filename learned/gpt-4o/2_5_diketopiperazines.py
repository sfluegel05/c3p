"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a piperazine-2,5-dione skeleton, which is a six-membered ring containing
    two nitrogen atoms and two carbonyl groups. Substitutions and stereochemistry are common.

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

    # Updated SMARTS pattern with stereochemistry considerations
    diketopiperazine_pattern = Chem.MolFromSmarts('O=C1NCC(=O)N1')
    if not mol.HasSubstructMatch(diketopiperazine_pattern):
        return False, "Does not contain piperazine-2,5-dione skeleton"

    # Use RingInfo to count the number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, f"Expected at least one ring; found {num_rings}."

    # Check for at least two oxygens typical for diketopiperazine scaffolds
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Too few oxygens for typical diketopiperazine: found {o_count}"

    # Consider typical sizes for 2,5-diketopiperazines, which can have diverse substitutions
    n_atoms = mol.GetNumAtoms()
    if n_atoms > 100:
        return False, "Molecule is too large to be a typical 2,5-diketopiperazine with complex substitutions"

    return True, "Contains piperazine-2,5-dione skeleton"