"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26976 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is derived from a tetraterpene (C40 backbone), possibly modified
    by rearrangement or removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetraterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude molecules containing heteroatoms not typical for tetraterpenoids (e.g., nitrogen, sulfur, halogens)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [1, 6, 8]:  # Allow only carbon, hydrogen, oxygen
            return False, f"Contains heteroatom ({atom.GetSymbol()}) not typical for tetraterpenoids"

    # Count number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 45:
        return False, f"Carbon count ({c_count}) not within typical range for tetraterpenoids"

    # Check for the presence of isoprene units
    # General isoprene unit pattern allowing for modifications
    isoprene_patterns = [
        Chem.MolFromSmarts("C(=C)C"),           # Isoprene unit
        Chem.MolFromSmarts("C(=C)CC"),          # Isoprene with rearrangement
        Chem.MolFromSmarts("C=C(C)C"),          # Isoprene alternative
        Chem.MolFromSmarts("C=C(C)CC"),         # Extended isoprene
        Chem.MolFromSmarts("C(=C)C(C)C"),       # Isoprene variant
    ]

    isoprene_match_count = 0
    for pattern in isoprene_patterns:
        matches = mol.GetSubstructMatches(pattern)
        isoprene_match_count += len(matches)

    if isoprene_match_count < 6:
        return False, f"Insufficient isoprene-like units found ({isoprene_match_count})"

    # Check for extended conjugation
    # Calculate the fraction of sp2 carbons
    sp2_carbons = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetHybridization() == Chem.HybridizationType.SP2:
            sp2_carbons += 1

    if sp2_carbons / c_count < 0.3:
        return False, "Insufficient sp2 hybridized carbons for extended conjugation"

    # Exclude molecules with unusual functional groups
    # Count number of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    # Allow some oxygen (e.g., hydroxyls, carbonyls), but limit excessive functionalization
    if o_count > 10:
        return False, f"Too many oxygen atoms ({o_count}) for a typical tetraterpenoid"

    # Calculate molecular weight
    mol_wt = Descriptors.ExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 1000:
        return False, f"Molecular weight ({mol_wt:.2f}) not typical for a tetraterpenoid"

    return True, "Molecule matches criteria for a tetraterpenoid (isoprenoid units, conjugation, typical functional groups)"