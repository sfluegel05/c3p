"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: chlorophyll

Definition: A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings. The rings can have various side chains which usually include a long phytol chain.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    A chlorophyll is characterized by a magnesium porphyrin core with a fifth ring and side chains, often including a long phytol chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium atom (atomic number 12)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found"

    # Check for porphyrin macrocycle (four pyrrole rings connected via methine bridges)
    porphyrin_smarts = 'C1=C[C@H]2C=C3N=C(C=C4N=CC(C=C1N2)=C34)'
    porphyrin_pattern = Chem.MolFromSmarts(porphyrin_smarts)
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin macrocycle found"

    # Check for fifth ring fused to porphyrin core (e.g., cyclopentanone ring)
    fifth_ring_smarts = 'C1=CC=CC2=CN=CC=C12'  # Simplified pattern for fused ring
    fifth_ring_pattern = Chem.MolFromSmarts(fifth_ring_smarts)
    if not mol.HasSubstructMatch(fifth_ring_pattern):
        return False, "No fifth ring fused to porphyrin core found"

    # Optionally, check for phytol side chain (long branched hydrocarbon chain)
    phytol_smarts = 'CC(C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C'
    phytol_pattern = Chem.MolFromSmarts(phytol_smarts)
    if not mol.HasSubstructMatch(phytol_pattern):
        return False, "No phytol side chain found"

    return True, "Contains magnesium porphyrin core with fifth ring and phytol side chain"

__metadata__ = {
    'chemical_class': {
        'name': 'chlorophyll',
        'definition': 'A family of magnesium porphyrins, defined by the presence of a fifth ring beyond the four pyrrole-like rings. The rings can have various side chains which usually include a long phytol chain.'
    }
}