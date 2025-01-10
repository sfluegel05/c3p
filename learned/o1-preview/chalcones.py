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
    A chalcone is a ketone that is 1,3-diphenylpropenone (benzylideneacetophenone), ArCH=CH(=O)Ar,
    and its derivatives formed by substitution.

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

    # Define chalcone core pattern: aromatic ring connected via alpha,beta-unsaturated ketone to another aromatic ring
    # Ensure the carbonyl carbon is not in a ring (R0) to exclude aurones
    chalcone_pattern = Chem.MolFromSmarts("[a]-[C;R0](=O)-[C;R0]=[C;R0]-[a]")
    # Define dihydrochalcone pattern: saturated version of chalcone
    dihydrochalcone_pattern = Chem.MolFromSmarts("[a]-[C;R0](=O)-[C;R0]-[C;R0]-[a]")

    # Check for the chalcone or dihydrochalcone core
    is_chalcone = mol.HasSubstructMatch(chalcone_pattern)
    is_dihydrochalcone = mol.HasSubstructMatch(dihydrochalcone_pattern)

    if not (is_chalcone or is_dihydrochalcone):
        return False, "Chalcone core not found"

    # Verify the number of aromatic rings (should be at least 2)
    ri = mol.GetRingInfo()
    num_aromatic_rings = 0
    for ring in ri.AtomRings():
        if all(mol.GetAtomWithIdx(atom_idx).GetIsAromatic() for atom_idx in ring):
            num_aromatic_rings +=1
    if num_aromatic_rings < 2:
        return False, f"Less than two aromatic rings detected ({num_aromatic_rings})"

    # Optionally, check that the core is acyclic (exclude ring systems involving the core)
    # Find all matches to the chalcone or dihydrochalcone pattern and ensure they are not in rings
    if is_chalcone:
        matches = mol.GetSubstructMatches(chalcone_pattern)
        pattern_name = "chalcone"
    else:
        matches = mol.GetSubstructMatches(dihydrochalcone_pattern)
        pattern_name = "dihydrochalcone"

    for match in matches:
        atoms_in_rings = any(mol.GetAtomWithIdx(idx).IsInRing() for idx in match)
        if atoms_in_rings:
            continue  # Skip if any atom in the match is in a ring
        else:
            return True, f"{pattern_name.capitalize()} core detected"

    return False, "Chalcone core not found or is part of a ring system"