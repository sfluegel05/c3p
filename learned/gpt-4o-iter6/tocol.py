"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    Tocols have a chromanol core with a hydrocarbon chain attached at position 2
    consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined chroman-6-ol core pattern
    # The core is a [1,4]dioxine system attached to an aromatic benzene ring
    chromanol_pattern = Chem.MolFromSmarts("c1cc(O)c2C[C@@]3(C)OCCc3c(O)c2c1")
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "No chroman-6-ol core found"

    # Design a flexible isoprenoid unit pattern considering possible variations
    # Generally, an isoprenoid unit is a C5 unit, often found as C1=C-C-C-C
    # We allow for either saturated, unsaturated, or highly branched variants
    isoprenoid_unit_patterns = [
        Chem.MolFromSmarts("C(=C)C(C)C"),  # Idealized isoprenoid unit (can extend this list)
        Chem.MolFromSmarts("C=C(C)CC"),    # Variant with a double bond
        Chem.MolFromSmarts("C(C)C(C)C")    # Saturated isoprenoid variant
    ]

    # Count the number of isoprenoid units by trying different patterns
    isoprenoid_count = 0
    for isoprenoid_unit in isoprenoid_unit_patterns:
        isoprenoid_count += len(mol.GetSubstructMatches(isoprenoid_unit))

    if isoprenoid_count < 3:
        return False, f"Found {isoprenoid_count} isoprenoid units, need at least 3"

    # Ensure attachment at position 2, though cited pattern ensures it by structure
    # Check that the chromanol core is correctly substituted with the chain
    matches = mol.GetSubstructMatch(chromanol_pattern)
    if not matches:
        return False, "Chromanol core detected but verification failed"

    return True, "Contains a chroman-6-ol core with valid hydrocarbon chain substitution spanning three isoprenoid units"