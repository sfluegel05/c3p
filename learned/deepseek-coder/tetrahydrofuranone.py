"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: CHEBI:XXXXX tetrahydrofuranone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is a tetrahydrofuran ring with an oxo group (C=O) attached to it.

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

    # Define the tetrahydrofuran ring pattern (C1CCOC1)
    tetrahydrofuran_pattern = Chem.MolFromSmarts("[C]1[C][C][O][C]1")
    if not mol.HasSubstructMatch(tetrahydrofuran_pattern):
        return False, "No tetrahydrofuran ring found"

    # Define the oxo group pattern (C=O)
    oxo_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if len(oxo_matches) == 0:
        return False, "No oxo group found"

    # Check if the oxo group is attached to the tetrahydrofuran ring
    # We need to ensure that the oxo group is directly connected to the ring
    ring_atoms = mol.GetRingInfo().AtomRings()
    if not ring_atoms:
        return False, "No ring found in the molecule"

    # Get the atoms in the tetrahydrofuran ring
    tetrahydrofuran_ring = None
    for ring in ring_atoms:
        if len(ring) == 5:  # Tetrahydrofuran is a 5-membered ring
            # Check if the ring contains an oxygen atom
            if any(mol.GetAtomWithIdx(atom).GetAtomicNum() == 8 for atom in ring):
                tetrahydrofuran_ring = ring
                break

    if tetrahydrofuran_ring is None:
        return False, "No tetrahydrofuran ring found"

    # Check if any oxo group is attached to the tetrahydrofuran ring
    for oxo_match in oxo_matches:
        oxo_carbon = oxo_match[0]
        if oxo_carbon in tetrahydrofuran_ring:
            return True, "Tetrahydrofuran ring with an oxo group attached"

    return False, "Oxo group not attached to the tetrahydrofuran ring"