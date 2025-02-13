"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: CHEBI:51584 tetrahydrofuranone
Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for tetrahydrofuran ring with an oxo substituent
    thf_oxo_pattern = Chem.MolFromSmarts("[O;R]1[C;R][C;R][C;R][C;R]1=O")
    thf_oxo_match = mol.HasSubstructMatch(thf_oxo_pattern)
    if not thf_oxo_match:
        return False, "No tetrahydrofuran ring with an oxo substituent found"

    # Check if the molecule contains only one ring
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 1:
        return False, "More than one ring present"

    # Check if the ring size is 5
    ring_atoms = ring_info.AtomRings()[0]
    if len(ring_atoms) != 5:
        return False, "Ring size is not 5"

    # Count heteroatoms in the ring
    hetero_count = sum(1 for atom in ring_atoms if atom.GetAtomicNum() != 6)
    if hetero_count != 1:
        return False, "Ring contains more than one heteroatom"

    return True, "Molecule is a tetrahydrofuranone"