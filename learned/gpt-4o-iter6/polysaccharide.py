"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem

def is_polysaccharide(smiles):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide consists of more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Define a more comprehensive pattern for monosaccharide units
    # Modified to potentially fit more common variations e.g. glucose unit
    cyclic_mono_pattern = Chem.MolFromSmarts("[C&R]1([O])[C&R]([C&R]([CO])[C&R]([O])[C&R]1O)")  # Pyranose-like structure
    acetamido_pattern = Chem.MolFromSmarts("NC(=O)C")  # Acetamido group common in GlcNAc

    # Iterate over potential substructures and count
    cycle_mono_count = len(mol.GetSubstructMatches(cyclic_mono_pattern))
    acetamido_count = len(mol.GetSubstructMatches(acetamido_pattern))

    # Consider both counts towards a theoretical mono count
    mono_count = cycle_mono_count + acetamido_count

    if mono_count < 10:
        return False, f"Contains {mono_count} probable monosaccharide rings/units, less than required 10 for polysaccharide"

    # Define pattern for glycosidic linkages
    glycosidic_link_pattern = Chem.MolFromSmarts("[C&R]O[C&R]")  # Ether-like linkage focusing on C-O-C

    # Count glycosidic linkages
    glyco_count = len(mol.GetSubstructMatches(glycosidic_link_pattern))
    if glyco_count < mono_count - 1:
        return False, f"Insufficient glycosidic linkages ({glyco_count}) for {mono_count} monosaccharide rings/units"

    # Final classification
    return True, "Contains sufficient monosaccharide moieties linked glycosidically to classify as a polysaccharide"