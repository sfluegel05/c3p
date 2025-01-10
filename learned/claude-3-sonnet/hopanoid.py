"""
Classifies: CHEBI:51963 hopanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    Hopanoids are triterpenoids based on a hopane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic hopane core pattern (simplified to match pentacyclic core)
    hopane_core = Chem.MolFromSmarts("[C]1~[C]~[C]~[C]2~[C]1~[C]~[C]3~[C]2~[C]~[C]4~[C]3~[C]~[C]5~[C]4~[C]~[C]~[C]5")
    if not mol.HasSubstructMatch(hopane_core):
        return False, "No hopane pentacyclic core found"

    # Count rings
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 5:
        return False, "Insufficient number of rings for hopanoid"

    # Count carbons and check range (typical hopanoids have 30-35 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25:  # Allow some flexibility but require substantial carbon skeleton
        return False, f"Too few carbons ({c_count}) for hopanoid structure"

    # Check for branching methyl groups (characteristic of hopanoids)
    methyl_pattern = Chem.MolFromSmarts("[CH3][C]")
    methyl_matches = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_matches < 4:  # Hopanoids typically have multiple methyl groups
        return False, f"Insufficient methyl groups ({methyl_matches}) for hopanoid"

    # Calculate molecular weight - hopanoids typically >400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight ({mol_wt:.1f}) too low for hopanoid"

    # Check for characteristic gem-dimethyl group often present in hopanoids
    gem_dimethyl = Chem.MolFromSmarts("[C]([CH3])([CH3])")
    if not mol.HasSubstructMatch(gem_dimethyl):
        return False, "Missing characteristic gem-dimethyl group"

    # Count sp3 carbons (hopanoids are highly saturated)
    sp3_carbons = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[CX4]")))
    if sp3_carbons < 20:
        return False, f"Insufficient sp3 carbons ({sp3_carbons}) for hopanoid structure"

    # Check for reasonable number of rings and their connectivity
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    if ring_count < 5 or ring_count > 7:  # Allow for additional small rings in derivatives
        return False, f"Unusual number of rings ({ring_count}) for hopanoid"

    # Success - molecule appears to be a hopanoid
    return True, "Contains hopane skeleton with characteristic features: pentacyclic core, " \
                 "appropriate methylation pattern, and molecular weight"