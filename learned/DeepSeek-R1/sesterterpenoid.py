"""
Classifies: CHEBI:26660 sesterterpenoid
"""
"""
Classifies: CHEBI:36509 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on SMILES string.
    Sesterterpenoids are derived from a C25 sesterterpene skeleton, possibly modified.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES"

    # Basic molecular formula check (C25 base with possible modifications)
    formula = Descriptors.CalcMolFormula(mol)
    c_count = sum(1 for c in formula if c == 'C')
    if not (20 <= c_count <= 30):  # Allow significant modification
        return False, f"C count {c_count} outside sesterterpenoid range"

    # Common terpenoid structural features
    terp_pattern = Chem.MolFromSmarts("""
        ([#6]-[#6](-[#6])-[#6](-[#6])-[#6] |  # Branched chain pattern
        C/C=C/\\C |                            # Isoprene-like double bonds
        [CH2]C(C)(C)                           # Terpenoid methyl branches
    """)
    if not mol.HasSubstructMatch(terp_pattern):
        return False, "Lacks terpenoid structural motifs"

    # Check for head-to-tail isoprene linkages (approximate)
    isoprene_linkage = Chem.MolFromSmarts("C/C=C/C(/C)=C/C")
    matches = len(mol.GetSubstructMatches(isoprene_linkage))
    if matches < 3:  # Expect multiple isoprene-like units
        return False, f"Only {matches} isoprene-like linkages found"

    # Look for cyclization patterns common in terpenoids
    ring_systems = Chem.GetSSSR(mol)
    if ring_systems < 1:  # Most sesterterpenoids have some rings
        return False, "No ring systems detected"

    # Check for oxidation patterns (common in terpenoids)
    oxygen_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if oxygen_count < 2:  # Most have some oxygenated groups
        return False, "Insufficient oxygen functional groups"

    return True, "Contains terpenoid features with C25-derived structure"