"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: CHEBI:36275 lactol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the lactol pattern: a ring with an oxygen and a hydroxyl group
    lactol_pattern = Chem.MolFromSmarts("[OH]-C1-O-C-C-C1")
    if not mol.HasSubstructMatch(lactol_pattern):
        return False, "No lactol pattern found (cyclic hemiacetal structure)"

    # Check for the presence of a carbonyl group (aldehyde or ketone) in the molecule
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found (aldehyde or ketone required for lactol formation)"

    # Verify that the molecule is cyclic
    if not mol.GetRingInfo().NumRings():
        return False, "Molecule is not cyclic (lactols are cyclic hemiacetals)"

    # Check for the presence of a hydroxyl group attached to the ring
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]-C1-O-C-C-C1")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group attached to the ring"

    # Check for the presence of an oxygen in the ring
    ring_oxygen_pattern = Chem.MolFromSmarts("C1-O-C-C-C1")
    if not mol.HasSubstructMatch(ring_oxygen_pattern):
        return False, "No oxygen in the ring (required for lactol structure)"

    return True, "Contains a cyclic hemiacetal structure with a hydroxyl group and an oxygen in the ring"

# Example usage:
# smiles = "C1C(C(C(O1)O)O)O"  # 2,3,4-trihydroxytetrahydrofuran
# print(is_lactol(smiles))