"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: CHEBI:3106 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    An alkaloid is a naturally occurring, basic nitrogen compound, usually heterocyclic,
    found in plants, bacteria, fungi, and animals. It excludes amino acids, peptides,
    proteins, nucleotides, nucleic acids, amino sugars, and antibiotics. Compounds with
    exocyclic nitrogen are usually classified as amines rather than alkaloids.

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

    # Check for basic nitrogen
    basic_nitrogen_smarts = "[N+;!X3;!X4]"
    basic_nitrogen_pattern = Chem.MolFromSmarts(basic_nitrogen_smarts)
    if not mol.HasSubstructMatch(basic_nitrogen_pattern):
        return False, "No basic nitrogen found"

    # Check for heterocyclic rings
    heterocyclic_ring_smarts = "[r;!R2]"  # Rings with at least one non-carbon atom
    heterocyclic_ring_pattern = Chem.MolFromSmarts(heterocyclic_ring_smarts)
    if not mol.HasSubstructMatch(heterocyclic_ring_pattern):
        return False, "No heterocyclic rings found"

    # Exclude common non-alkaloid compounds
    excluded_smarts = (
        "[NX3;H2,H1]",  # Amino acids
        "[NX4;H3,H2,H1]",  # Amines
        "[C&r3,r4,r5,r6]",  # Peptides, proteins
        "[n;r3,r4,r5,r6]",  # Nucleotides, nucleic acids
        "[O&r3,r4,r5,r6]",  # Sugars
        "[C&r3,r4,r5,r6][N&r3,r4,r5,r6]",  # Peptides, proteins
    )
    for smarts in excluded_smarts:
        pattern = Chem.MolFromSmarts(smarts)
        if mol.HasSubstructMatch(pattern):
            return False, "Excluded compound found"

    # Count aromatic rings and rings with basic nitrogen
    aromatic_rings = mol.GetRingInfo().AtomRings()
    basic_nitrogen_rings = []
    for ring in aromatic_rings:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        if any(atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1 for atom in ring_atoms):
            basic_nitrogen_rings.append(ring)

    # Alkaloids typically have 1-3 aromatic rings, at least one with a basic nitrogen
    if not (1 <= len(aromatic_rings) <= 3) or len(basic_nitrogen_rings) == 0:
        return False, "Aromatic ring count or basic nitrogen ring count outside expected range"

    # Check molecular weight - alkaloids typically <1000 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 1000:
        return False, "Molecular weight too high for alkaloid"

    # Passed all checks, classify as alkaloid
    return True, "Contains basic nitrogen, heterocyclic rings, and meets other structural requirements"