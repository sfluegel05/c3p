"""
Classifies: CHEBI:18154 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Define a general pattern for cyclic monosaccharide units (e.g., hexose)
    cyclic_mono_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C1O")  # Simplified pyranose form
    acetamido_pattern = Chem.MolFromSmarts("N(C(*=O)C)C")  # Acetamido group indicating sugar derivatives like GlcNAc
    
    # Count cyclic monosaccharide units
    cycle_mono_count = len(mol.GetSubstructMatches(cyclic_mono_pattern))
    acetamido_count = len(mol.GetSubstructMatches(acetamido_pattern))

    # Consider both counts towards a theoretical mono count
    mono_count = cycle_mono_count + acetamido_count

    if mono_count < 10:
        return False, f"Contains {mono_count} probable monosaccharide rings/units, less than required 10 for polysaccharide"

    # Glycosidic linkage pattern (C-O-C bond excluding terminal connections)
    glycosidic_link_pattern = Chem.MolFromSmarts("[!#1]O[#6&R]([!#1])")  # Ether-like O-C linkage

    # Count glycosidic linkages
    glyco_count = len(mol.GetSubstructMatches(glycosidic_link_pattern))
    if glyco_count < mono_count - 1:
        return False, f"Insufficient glycosidic linkages ({glyco_count}) for {mono_count} monosaccharide rings/units"

    # Final classification
    return True, "Contains sufficient monosaccharide moieties linked glycosidically to classify as a polysaccharide"