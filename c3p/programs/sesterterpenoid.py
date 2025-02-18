"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    A sesterterpenoid is derived from a sesterterpene, showcasing a modified C25 backbone, 
    rearranged or involving additional ring structures and diverse functional groups.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a sesterterpenoid, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Widely allow for structural variability in sesterterpenoids
    if c_count < 20 or c_count > 40:
        return False, f"Carbon count of {c_count} is outside the expanded range typical for sesterterpenoids"

    # Check for rings; sesterterpenoids frequently have multiple rings
    ring_count = mol.GetRingInfo().NumRings()
    if ring_count < 1:
        return False, "Contains less than 1 ring, uncommon for sesterterpenoids"

    # Check for modified isoprene-like structures or large frameworks
    modified_structure_count = 0
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)C")
    for match in mol.GetSubstructMatches(isoprene_pattern):
        modified_structure_count += 1

    if modified_structure_count < 2:
        return False, "Lacks enough recognized isoprene-like structural units for sesterterpenoids"

    # Look for specific terpenoid-like motifs: cyclohexane or cyclopentane units
    cyclic_motif_pattern = Chem.MolFromSmarts("C1CCCCC1 | C1CCCC1")
    if len(mol.GetSubstructMatches(cyclic_motif_pattern)) < 1:
        return False, "Missing key cyclic motifs common in terpenoids"

    # Functional group diversity; check for common sesterterpenoid traits
    oxygen_pattern = Chem.MolFromSmarts("[OH]")
    oxygen_matches = len(mol.GetSubstructMatches(oxygen_pattern))
    
    if oxygen_matches == 0:
        return False, "No hydroxyl groups which are often present in sesterterpenoid moieties"

    # If all checks are satisfied, classify as sesterterpenoid
    return True, "Molecule exhibits a modified terpenoid backbone with requisite structural diversity typical of a sesterterpenoid."