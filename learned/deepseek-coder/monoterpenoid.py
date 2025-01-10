"""
Classifies: CHEBI:25409 monoterpenoid
"""
"""
Classifies: CHEBI:35189 monoterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid based on its SMILES string.
    A monoterpenoid is derived from a monoterpene (C10 skeleton) and may have rearrangements or modifications.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Monoterpenoids typically have 10 carbons, but modifications may reduce this
    if c_count < 8 or c_count > 12:
        return False, f"Carbon count ({c_count}) is outside the expected range for a monoterpenoid"

    # Check for isoprene-like patterns (C5 units)
    isoprene_pattern = Chem.MolFromSmarts("[CH3]-[CH2]-[CH]=[CH2]")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene-like pattern found"

    # Check for cyclic or acyclic structures typical of monoterpenoids
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings > 2:
        return False, f"Too many rings ({num_rings}) for a monoterpenoid"

    # Check for functional groups common in monoterpenoids (e.g., alcohols, ketones, esters)
    functional_groups = ["[OH]", "[C=O]", "[O][C=O]"]
    has_functional_group = any(mol.HasSubstructMatch(Chem.MolFromSmarts(fg)) for fg in functional_groups)
    if not has_functional_group:
        return False, "No typical monoterpenoid functional groups found"

    # Check molecular weight (monoterpenoids typically have MW between 130 and 200)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 130 or mol_wt > 200:
        return False, f"Molecular weight ({mol_wt:.2f}) is outside the expected range for a monoterpenoid"

    return True, "Contains a C10 skeleton or modified form with typical monoterpenoid features"