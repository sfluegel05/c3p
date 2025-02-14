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
    which may have rearranged or modified skeletons, and can have various functional groups.

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

    # Check for presence of elements typical in diterpenoids
    allowed_elements = {6, 1, 7, 8, 15, 16}  # C, H, N, O, P, S
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_elements:
            return False, f"Contains element {atom.GetSymbol()} not typical in diterpenoids"

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Diterpenoids generally have around 20 carbons but can vary due to modifications
    if c_count < 15 or c_count > 45:
        return False, f"Carbon count ({c_count}) not typical for diterpenoids (15-45 carbons)"

    # Calculate Double Bond Equivalents (DBE)
    dbe = rdMolDescriptors.CalcNumRings(mol) + mol.GetNumBonds() - mol.GetNumAtoms() + 1

    # Diterpenoids often have multiple rings
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 2:
        return False, f"Too few rings ({num_rings}) for typical diterpenoids"

    # Check for common diterpenoid skeletons using SMARTS patterns
    diterpene_scaffolds = [
        # Abietane skeleton
        Chem.MolFromSmarts('C1(CC[C@@]2([C@@H]1CC=C2)C)C'),
        # Labdane skeleton
        Chem.MolFromSmarts('C1CC2CCC1(C)C(CCC2C)C'),
        # Clerodane skeleton
        Chem.MolFromSmarts('C1CC2CCC1(C)C=CC2'),
        # Additional diterpene skeletons can be added
    ]

    scaffold_matches = False
    for scaffold in diterpene_scaffolds:
        if mol.HasSubstructMatch(scaffold):
            scaffold_matches = True
            break

    if not scaffold_matches:
        return False, "Does not contain common diterpenoid skeletons"

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
    return True, "Molecule has features consistent with diterpenoids (common skeleton, ring count)"