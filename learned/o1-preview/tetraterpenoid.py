"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26964 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    A tetraterpenoid is a terpenoid derived from a tetraterpene (C40 skeleton).
    It typically has eight isoprene units connected head-to-tail, forming a long
    conjugated polyene chain, possibly with rings and functional groups.

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

    # Count number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35 or c_count > 45:
        return False, f"Carbon count ({c_count}) not in range for tetraterpenoid (35-45 carbons)"
    
    # Check for long conjugated polyene chain (10 or more alternating double bonds)
    # Create a pattern for conjugated double bonds
    polyene_pattern = Chem.MolFromSmarts("C(=C)C=C" * 5)  # Pattern for 10 conjugated double bonds
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No long conjugated polyene chain found (need at least 10 alternating double bonds)"
    
    # Check for isoprene units (C5 units)
    # Isoprene unit pattern: C=C-C-C=C
    isoprene_pattern = Chem.MolFromSmarts("C=C-C-C=C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 6:
        return False, f"Found {len(isoprene_matches)} isoprene units, need at least 6"
    
    # Check for presence of methyl substituents on the polyene chain
    methyl_pattern = Chem.MolFromSmarts("[CH3]-[C]=[C]")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) < 4:
        return False, f"Found {len(methyl_matches)} methyl substituents on polyene chain, need at least 4"
    
    # Check for terpenoid functional groups (oxygen-containing groups)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found; may not be a terpenoid"

    # Optional: Check for rings (some tetraterpenoids have rings)
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count > 6:
        return False, f"Too many rings ({ring_count}); atypical for tetraterpenoids"

    return True, "Molecule meets criteria for tetraterpenoid (C40-derived terpenoid with characteristic features)"