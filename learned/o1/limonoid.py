"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a triterpenoid that is highly oxygenated and has a prototypical structure
    containing or derived from a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1: Check for the steroid nucleus with 4 fused rings (6-6-6-5) and specific methyl groups
    # SMARTS pattern for steroid nucleus with 4,4,8-trimethyl substitutions
    steroid_pattern = Chem.MolFromSmarts("""
        [#6]1([#6])[#6][#6]2[#6]([#6])[#6][#6]3[C;H]([C;H](C)[C;H](C)[C;H]12)[#6][#6]4[C;H]3[C;H](C)[C;H](C)[C;H]4C
    """)
    if steroid_pattern is None:
        return False, "Invalid SMARTS pattern for steroid nucleus"

    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid nucleus with 4,4,8-trimethyl substitutions found"

    # Step 2: Check for furan ring attached at position 17
    # SMARTS pattern for furan ring connected to steroid nucleus
    furan_attachment_pattern = Chem.MolFromSmarts("""
        [#6]1([#6])[#6][#6]2[#6]([#6])[#6][#6]3[C;H]([C;H](C)[C;H](C)[C;H]12)[#6][#6]4[C;H]3[C;H](C)[C;H](C)[C;H]4C[$([#6][#6]c1ccoc1)]
    """)
    if furan_attachment_pattern is None:
        return False, "Invalid SMARTS pattern for furan attachment"

    if not mol.HasSubstructMatch(furan_attachment_pattern):
        return False, "No furan ring attached to steroid nucleus at position 17"

    # Step 3: Check if molecule is highly oxygenated (e.g., at least 6 oxygen atoms)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, f"Not highly oxygenated, contains only {o_count} oxygen atoms"

    # Step 4: Look for common functional groups in limonoids (e.g., lactones, acetates)
    lactone_pattern = Chem.MolFromSmarts('C(=O)O[C;H]')
    acetate_pattern = Chem.MolFromSmarts('C(=O)O[C;H]C')
    has_lactone = mol.HasSubstructMatch(lactone_pattern)
    has_acetate = mol.HasSubstructMatch(acetate_pattern)
    if not (has_lactone or has_acetate):
        return False, "No common limonoid functional groups (lactones or acetates) found"

    return True, "Molecule is a limonoid (contains 4,4,8-trimethyl-17-furanylsteroid skeleton and is highly oxygenated)"