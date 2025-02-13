"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: chalcones
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone that is 1,3-diphenylpropenone (benzylideneacetophenone),
    ArC(=O)CH=CHAr, and its derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define chalcone core pattern: aromatic ring connected via α,β-unsaturated ketone to another aromatic ring
    chalcone_pattern = Chem.MolFromSmarts("[a][C](=O)[C]=[C][a]")
    # Define dihydrochalcone pattern: saturated version of chalcone
    dihydrochalcone_pattern = Chem.MolFromSmarts("[a][C](=O)[C][C][a]")

    # Check for the chalcone or dihydrochalcone core
    is_chalcone = mol.HasSubstructMatch(chalcone_pattern)
    is_dihydrochalcone = mol.HasSubstructMatch(dihydrochalcone_pattern)

    if not (is_chalcone or is_dihydrochalcone):
        return False, "Chalcone core not found"

    # Verify the number of aromatic rings (should be at least 2)
    aromatic_rings = []
    ri = mol.GetRingInfo()
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring):
            aromatic_rings.append(ring)
    num_aromatic_rings = len(aromatic_rings)
    if num_aromatic_rings < 2:
        return False, f"Less than two aromatic rings detected ({num_aromatic_rings})"

    # Find all matches to the chalcone or dihydrochalcone pattern
    if is_chalcone:
        matches = mol.GetSubstructMatches(chalcone_pattern)
        pattern_name = "chalcone"
    else:
        matches = mol.GetSubstructMatches(dihydrochalcone_pattern)
        pattern_name = "dihydrochalcone"

    if matches:
        return True, f"{pattern_name.capitalize()} core detected"
    else:
        return False, "Chalcone core not found"