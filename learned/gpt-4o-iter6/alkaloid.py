"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    Alkaloids are naturally occurring, basic nitrogen compounds (often heterocyclic)
    found in living organisms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitrogen
    if not any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms()):
        return False, "Molecule does not contain nitrogen"

    # Check for nitrogen in ring structures
    ring_info = mol.GetRingInfo()
    has_nitrogen_in_ring = False
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
            has_nitrogen_in_ring = True
            break

    if not has_nitrogen_in_ring:
        return False, "Nitrogen is not part of a heterocyclic structure"

    # Detect extended set of broad nitrogen-containing structures
    broad_patterns = [
        "[nR]",        # Generic heteroaromatic rings
        "[NX3;R]",     # Tertiary amine in rings
        "[NX2;R]"      # Secondary amine in rings
    ]
    
    has_broad_heterocycle = any(mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)) for pattern in broad_patterns)
    if not has_broad_heterocycle:
        return False, "Structure doesn't match broad heterocyclic alkaloid pattern"

    # Allow flexibility in peptide-like structures
    peptide_like_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    peptide_matches = mol.GetSubstructMatches(peptide_like_pattern)
    if len(peptide_matches) > 4:
        # Some flexibility on amide patterns to allow diverse alkaloids
        return False, "Appears strongly peptide-like rather than typical alkaloid"

    return True, "Contains nitrogen in a heterocyclic structure typical of alkaloids"