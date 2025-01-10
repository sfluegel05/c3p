"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings found"

    # Find largest ring size
    max_ring_size = max(len(ring) for ring in ring_info.AtomRings())
    if max_ring_size < 12:
        return False, f"Largest ring has {max_ring_size} members, need at least 12"

    # Look for lactone group (cyclic ester)
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[OX2][#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"

    # Check if ester is part of a ring
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ester_in_ring = any(any(atom_idx in ring for ring in ring_info.AtomRings())
                       for match in ester_matches for atom_idx in match)
    if not ester_in_ring:
        return False, "Ester group not part of a ring"

    # Check for polyketide-like structure (alternating carbonyl and hydroxyl groups)
    polyketide_pattern = Chem.MolFromSmarts("[#6][CX3](=O)[#6][OX2H]")
    if not mol.HasSubstructMatch(polyketide_pattern):
        return False, "No polyketide-like structure found"

    # Additional checks for typical macrolide features
    # Count oxygen atoms (macrolides typically have multiple oxygens)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Too few oxygen atoms for macrolide"

    # Check molecular weight (macrolides are typically >300 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for macrolide"

    return True, "Contains macrocyclic lactone ring with 12+ members and polyketide-derived structure"