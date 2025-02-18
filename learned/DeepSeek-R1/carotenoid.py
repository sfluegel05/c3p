"""
Classifies: CHEBI:23044 carotenoid
"""
"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are tetraterpenoids (C40-derived) with conjugated polyene chains, excluding retinoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check basic carbon skeleton requirements
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 35:  # Allow some modification from C40 base
        return False, f"Insufficient carbons ({c_count} < 35)"

    # Detect conjugated polyene chain (at least 8 alternating double bonds)
    polyene_pattern = Chem.MolFromSmarts("*=*-#*=*-#*=*-#*=*-#*=*-#*=*-#*=*-#*=*")
    if not polyene_pattern:
        return None, None  # Fallback if pattern fails to compile
    if not mol.HasSubstructMatch(polyene_pattern):
        return False, "No extended conjugated polyene chain"

    # Check for isoprenoid methyl groups (at least 4 methyls near double bonds)
    methyl_pattern = Chem.MolFromSmarts("[CH3]-[C]=[C]")
    if not methyl_pattern:
        return None, None
    methyl_matches = len(mol.GetSubstructMatches(methyl_pattern))
    if methyl_matches < 4:
        return False, f"Only {methyl_matches} isoprenoid methyl groups"

    # Exclude retinoids using beta-ionone ring check
    beta_ionone = Chem.MolFromSmarts("[C]1([CH3])=C([CH3])[CH2][CH2][C]([CH3])([CH3])[CH2]1")
    if beta_ionone and mol.HasSubstructMatch(beta_ionone):
        if c_count < 40:  # Retinoids have shorter chains
            return False, "Beta-ionone ring with short chain (retinoid)"

    # Molecular weight range (400-1000 Da typical for carotenoids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} Da outside typical range"

    # Check oxygen count (carotenoids can have 0-6 oxygens typically)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 6:
        return False, f"Too many oxygens ({o_count}) for carotenoid"

    return True, "Extended conjugated polyene with isoprenoid features"