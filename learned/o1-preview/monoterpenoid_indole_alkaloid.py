"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
"""
Classifies: monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid contains an indole moiety derived from L-tryptophan and
    a monoterpene-derived moiety (C10 unit), usually linked together.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for indole moiety (both protonated and unprotonated nitrogen)
    indole_pattern = Chem.MolFromSmarts('c1ccc2c(c1)[nH]cc2')  # Indole with [nH]
    indole_pattern_unprot = Chem.MolFromSmarts('c1ccc2c(c1)ncc2')  # Indole with nitrogen

    # Check for indole moiety
    if not mol.HasSubstructMatch(indole_pattern) and not mol.HasSubstructMatch(indole_pattern_unprot):
        return False, "No indole moiety found"

    # Check for monoterpene-derived moiety (approximate)
    # Monoterpenes are C10 units built from isoprene units (C5H8)
    # Attempt to find two isoprene units connected to indole moiety

    # Count total number of carbons
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 20:  # Adjusted threshold to reduce false positives
        return False, f"Too few carbons ({num_carbons}) to be a monoterpenoid indole alkaloid"

    # Check for minimum number of rings (monoterpenoid indole alkaloids are typically polycyclic)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Too few rings ({num_rings}) to be a monoterpenoid indole alkaloid"

    # Attempt to find isoprene units connected to indole
    # Isoprene unit SMARTS pattern: C=C-C=C
    isoprene_pattern = Chem.MolFromSmarts('C=C-C=C')
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, "Less than two isoprene units found"

    # Check for connectivity between indole and monoterpene moieties
    # Attempt to find path between indole nitrogen and isoprene units
    indole_matches = mol.GetSubstructMatches(indole_pattern) + mol.GetSubstructMatches(indole_pattern_unprot)
    indole_atoms = [atom_idx for match in indole_matches for atom_idx in match]
    isoprene_atoms = [atom_idx for match in isoprene_matches for atom_idx in match]

    # Create a set of all atoms in indole and isoprene units
    indole_atom_set = set(indole_atoms)
    isoprene_atom_set = set(isoprene_atoms)

    # Check if there is a path connecting indole and isoprene atoms
    paths_found = False
    for indole_atom_idx in indole_atom_set:
        for isoprene_atom_idx in isoprene_atom_set:
            path = Chem.rdmolops.GetShortestPath(mol, indole_atom_idx, isoprene_atom_idx)
            if path and len(path) <= 5:  # Adjust path length as needed
                paths_found = True
                break
        if paths_found:
            break

    if not paths_found:
        return False, "No linkage between indole and monoterpene units found"

    # Passed all checks
    return True, "Contains indole moiety connected to monoterpene-derived units"

__metadata__ = {
    'chemical_class': {
        'name': 'monoterpenoid indole alkaloid',
        'definition': 'A terpenoid indole alkaloid which is biosynthesised from L-tryptophan and diisoprenoid (usually secologanin) building blocks.'
    },
    'config': {
        # Configuration parameters if any
    }
}