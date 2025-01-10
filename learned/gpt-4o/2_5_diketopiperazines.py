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

    # More flexible pattern allowing for substitutions on the diketopiperazine backbone
    diketopiperazine_pattern = Chem.MolFromSmarts('C1C(=O)NCC(=O)N1')
    if not mol.HasSubstructMatch(diketopiperazine_pattern):
        return False, "Does not contain piperazine-2,5-dione skeleton"

    # Verify that the skeleton is part of a larger chemistry space typical for these compounds
    n_rings = Chem.GetSSSR(mol)
    if n_rings < 1:
        return False, f"Expected ring systems; found {n_rings}."

    # Check additional properties to reduce false positives; diketopiperazines should have at least two oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Too few oxygens for typical diketopiperazine: found {o_count}"

    # Consider check for expected size or presence of typical functional groups
    if len(mol.GetAtoms()) > 50:
        return False, "Molecule is too large to be typical diketopiperazine without complex substitutions"

    return True, "Contains piperazine-2,5-dione skeleton"