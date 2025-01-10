"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a highly oxygenated triterpenoid with a significant structure containing or derived from a 4,4,8-trimethyl-17-furanylsteroid skeleton.
    
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

    # Check for carbon structure (considering broad variation typical for limonoids)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 23 or c_count > 40:
        return False, f"Number of carbon atoms {c_count} outside expected range for limonoid"

    # Ensure high oxygen content
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, f"Insufficient oxygen content for a highly oxygenated limonoid: {o_count}"

    # Look for key functional groups common in limonoids
    oxygen_patterns = [Chem.MolFromSmarts("C=O"), Chem.MolFromSmarts("C(=O)O"), Chem.MolFromSmarts("[OX2H]"), Chem.MolFromSmarts("O1C=CO1")]
    if not any(mol.HasSubstructMatch(pat) for pat in oxygen_patterns):
        return False, "Missing characteristic oxygenated functionalities (e.g., ketone, lactone, hydroxyl, furan)"

    # Check specific limonoid structural features - notably for a furan ring
    furan_ring_pattern = Chem.MolFromSmarts("o1ccc1")
    if not mol.HasSubstructMatch(furan_ring_pattern):
        return False, "Lack of furan ring structure which is common in limonoids"

    trimethyl_pattern = Chem.MolFromSmarts("C(C)(C)C")
    if not mol.HasSubstructMatch(trimethyl_pattern):
        return False, "Lack of characteristic trimethyl group for limonoid structure"

    return True, "SMILES string corresponds to a limonoid structure"