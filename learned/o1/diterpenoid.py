"""
Classifies: CHEBI:23849 diterpenoid
"""
"""
Classifies: CHEBI:35196 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid is derived from a diterpene (originally composed of four isoprene units),
    but may have rearranged or modified skeletons, often missing or adding skeletal atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate molecular formula
    formula = rdMolDescriptors.CalcMolFormula(mol)
    # Extract counts of elements
    from collections import Counter
    import re

    # Use regular expressions to extract element counts from the formula
    elements = re.findall('([A-Z][a-z]*)(\d*)', formula)
    element_counts = Counter()
    for elem, count in elements:
        count = int(count) if count else 1
        element_counts[elem] += count

    c_count = element_counts.get('C', 0)
    h_count = element_counts.get('H', 0)
    n_count = element_counts.get('N', 0)
    o_count = element_counts.get('O', 0)
    s_count = element_counts.get('S', 0)
    p_count = element_counts.get('P', 0)
    other_elements = set(element_counts.keys()) - {'C', 'H', 'O', 'N', 'S', 'P'}

    # Diterpenoids are generally composed of carbon, hydrogen, oxygen, nitrogen, sulfur, phosphorus
    if other_elements:
        return False, f"Contains elements not typical in diterpenoids: {', '.join(other_elements)}"

    # Diterpenoids generally have carbon counts around 20 (from four isoprene units)
    # But due to modifications, the carbon count can vary
    if c_count < 15 or c_count > 40:
        return False, f"Carbon count ({c_count}) not typical for diterpenoids (15-40 carbons)"

    # Calculate Double Bond Equivalents (DBE)
    dbe = (2 * c_count + 2 + n_count - h_count - halogen_count(mol)) / 2
    # Diterpenoids are generally polycyclic compounds with several rings
    if dbe < 4:
        return False, f"DBE ({dbe}) too low for diterpenoids"

    # Check for atypical functional groups
    atypical_groups = [
        Chem.MolFromSmarts('[N+](=O)[O-]'),  # Nitro group
        Chem.MolFromSmarts('S(=O)(=O)[O-]'),  # Sulfate group
        Chem.MolFromSmarts('P(=O)(O)(O)O'),  # Phosphate group
        Chem.MolFromSmarts('[Si]'),          # Silicon containing groups
        Chem.MolFromSmarts('[B]'),           # Boron containing groups
        Chem.MolFromSmarts('[F,Cl,Br,I]'),   # Halogens
    ]
    for group in atypical_groups:
        if mol.HasSubstructMatch(group):
            return False, "Contains atypical functional groups for diterpenoids"

    # If the molecule passes all checks, classify as diterpenoid
    return True, "Molecule has features consistent with diterpenoids (element composition, carbon count, DBE)"

def halogen_count(mol):
    """
    Counts the number of halogen atoms in the molecule.

    Args:
        mol (rdkit.Chem.Mol): Molecule object

    Returns:
        int: Number of halogen atoms
    """
    halogens = [9, 17, 35, 53, 85]  # Atomic numbers for F, Cl, Br, I, At
    return sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in halogens)