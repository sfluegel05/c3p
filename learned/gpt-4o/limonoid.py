"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a highly oxygenated triterpenoid containing a trimethyl-17-furanylsteroid skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for triterpenoid skeleton: 30 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 30:
        return False, f"Incorrect number of carbon atoms for a triterpenoid: {c_count}"

    # Check for presence of furanyl group
    furanyl_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furanyl_pattern):
        return False, "No furanyl group found"

    # Check for high oxygen content (limonoids are highly oxygenated)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, f"Oxygen count too low for limonoid: {o_count}"

    # Additional check for characteristic functionalities
    # Check for lactone group
    lactone_pattern = Chem.MolFromSmarts("C(=O)O")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone group found"

    # Check for presence of methyl groups which are common in trimethyl structure
    methyl_group_pattern = Chem.MolFromSmarts("C(C)(C)C")
    methyl_matches = [match for match in mol.GetSubstructMatches(methyl_group_pattern) if len(match) == 4]
    if len(methyl_matches) < 3:
        return False, f"Insufficient methyl groups for trimethyl structure: {len(methyl_matches)} groups"

    return True, "SMILES string corresponds to a limonoid structure"