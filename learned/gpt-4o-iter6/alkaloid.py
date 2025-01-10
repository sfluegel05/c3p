"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is characterized by having a basic nitrogen atom within a heterocyclic ring,
    and typically does not resemble a peptide or related compound.

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

    # Check for the presence of nitrogen
    has_nitrogen = any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    if not has_nitrogen:
        return False, "Molecule does not contain nitrogen"

    # Check for heterocyclic nitrogen
    ring_info = mol.GetRingInfo()
    has_nitrogen_in_ring = False
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
            has_nitrogen_in_ring = True
            break

    if not has_nitrogen_in_ring:
        return False, "Nitrogen is not part of a heterocyclic structure"

    # Detect broad heteroaromatic nitrogen including indoles/pyrrolidines/quinolines etc.
    broad_heteroaromatic_pattern = Chem.MolFromSmarts("[nR]")
    if not mol.HasSubstructMatch(broad_heteroaromatic_pattern):
        return False, "Structure doesn't match broad heterocyclic alkaloid pattern"

    # Avoid classifying peptide-like structures while allowing amide diversity
    peptide_like_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    peptide_matches = mol.GetSubstructMatches(peptide_like_pattern)
    if len(peptide_matches) > 3:
        return False, "Appears more like a peptide or related compound"

    return True, "Contains nitrogen in a heterocyclic structure typical of alkaloids"