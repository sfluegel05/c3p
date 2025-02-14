"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a biomacromolecule consisting of more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove explicit hydrogens for simplification
    mol = Chem.RemoveHs(mol)

    # Define a SMARTS pattern for monosaccharide units (pyranose and furanose rings with exocyclic oxygen)
    # This pattern looks for 5 or 6 membered rings containing oxygen and exocyclic oxygen (glycosidic linkage)
    monosaccharide_pattern = Chem.MolFromSmarts(
        "[C;R][C;R][C;R][C;R][O;R][C;R]"  # Pyranose ring
        "~[$([OX2H]),$([OX1]-C)]"         # Exocyclic oxygen (glycosidic linkage)
    )

    if monosaccharide_pattern is None:
        return False, "Error in SMARTS pattern definition"

    # Find unique matches for monosaccharide units
    mono_matches = mol.GetSubstructMatches(monosaccharide_pattern, useChirality=True)
    # Convert to sets of atom indices to ensure uniqueness
    mono_rings = [set(match[:6]) for match in mono_matches]  # First 6 atoms are the ring

    # Remove duplicate rings
    unique_rings = []
    for ring in mono_rings:
        if ring not in unique_rings:
            unique_rings.append(ring)

    total_units = len(unique_rings)

    if total_units > 10:
        return True, f"Contains {total_units} monosaccharide units linked glycosidically"
    else:
        return False, f"Contains only {total_units} monosaccharide units, need more than 10 for a polysaccharide"