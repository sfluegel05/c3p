"""
Classifies: CHEBI:26935 tetraterpenoid
"""
"""
Classifies: CHEBI:26764 tetraterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetraterpenoid(smiles: str):
    """
    Determines if a molecule is a tetraterpenoid based on its SMILES string.
    Tetraterpenoids are derived from tetraterpenes (C40 skeleton) and may have 
    modifications like rearrangements or removal of methyl groups.

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

    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Count methyl groups
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    methyl_count = len(mol.GetSubstructMatches(methyl_pattern))
    
    # Count double bonds
    double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    
    # Count rings
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    
    # Count oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    # Check for characteristic carotenoid backbone patterns
    carotenoid_patterns = [
        "C=CC=CC=CC=CC=C",  # Long conjugated chain
        "C=CC=CC=CC=CC=CC=C",  # Extended conjugation
        "C(C)(C)C=CC=CC=C",  # Typical carotenoid fragment
        "C1C(C)=CCCC1(C)C",  # Cyclic end group
        "CC(C)=CCCC(C)=CC=C" # Acyclic end with methyls
    ]
    
    backbone_matches = 0
    for pattern in carotenoid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            backbone_matches += 1

    # Check for nitrogen (rare in tetraterpenoids)
    if sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7) > 0:
        return False, "Contains nitrogen, unusual for tetraterpenoid"

    # Base molecule checks (for non-modified tetraterpenoids)
    if 38 <= c_count <= 42 and methyl_count >= 8 and double_bond_count >= 9 and backbone_matches >= 2:
        return True, f"Matches typical tetraterpenoid pattern with {c_count} carbons, {methyl_count} methyl groups, {double_bond_count} double bonds"

    # Modified tetraterpenoid checks
    is_modified = False
    reason = ""

    # Check for apo-carotenoids (shortened derivatives)
    if 20 <= c_count < 38 and methyl_count >= 4 and double_bond_count >= 5 and backbone_matches >= 1:
        if mol.HasSubstructMatch(Chem.MolFromSmarts("[CH]=O")) or \
           mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)O[H]")):
            is_modified = True
            reason = "Appears to be an apo-carotenoid (shortened tetraterpenoid)"

    # Check for glycosylated derivatives
    elif 42 < c_count <= 55 and o_count >= 5 and backbone_matches >= 2:
        if mol.HasSubstructMatch(Chem.MolFromSmarts("OC1C(O)C(O)C(O)C(O)C1")):  # Sugar pattern
            is_modified = True
            reason = "Appears to be a glycosylated tetraterpenoid"

    # Check for prephytoene-type compounds
    elif 38 <= c_count <= 42 and methyl_count >= 8 and \
         mol.HasSubstructMatch(Chem.MolFromSmarts("COP(=O)(O)OP(=O)(O)O")):
        is_modified = True
        reason = "Appears to be a prephytoene-type tetraterpenoid"

    # Check for other modified forms with characteristic features
    elif 20 <= c_count <= 55 and methyl_count >= 4 and backbone_matches >= 1:
        if double_bond_count >= 5 and (
            mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)C")) or  # Keto group
            mol.HasSubstructMatch(Chem.MolFromSmarts("C(O)C")) or   # Hydroxy group
            mol.HasSubstructMatch(Chem.MolFromSmarts("C1OC1")) or   # Epoxide
            mol.HasSubstructMatch(Chem.MolFromSmarts("C#C"))        # Triple bond
        ):
            is_modified = True
            reason = "Appears to be a modified tetraterpenoid with characteristic functional groups"

    if is_modified:
        return True, f"{reason} ({c_count} carbons, {methyl_count} methyl groups, {double_bond_count} double bonds)"

    return False, f"Does not match tetraterpenoid patterns (C:{c_count}, Me:{methyl_count}, DB:{double_bond_count})"