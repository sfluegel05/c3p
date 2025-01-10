"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of nitrogen
    has_nitrogen = any(atom.GetAtomicNum() == 7 for atom in mol.GetAtoms())
    if not has_nitrogen:
        return False, "Molecule does not contain nitrogen"

    # Check for heterocyclic structure
    # A basic check is to see if there is a ring that contains nitrogen
    ring_info = mol.GetRingInfo()
    has_nitrogen_in_ring = False
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
            has_nitrogen_in_ring = True
            break

    if not has_nitrogen_in_ring:
        return False, "Nitrogen is not part of a heterocyclic structure"

    # Check basic nature
    # This is tricky to identify directly from SMILES without additional computational methods.
    # Here, we assume if it has nitrogen in heterocyclic format, it might be basic.
    
    # Exclude amino acids, peptides, etc. (complex to do with SMILES alone)
    # We take a conservative approach and assume certain features might indicate these exclusions.
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num == 8:  # Oxygen-heavy might indicate non-alkaloid features
            return False, "Appears more like a peptide or related compound"

    return True, "Contains nitrogen in a heterocyclic structure typical of alkaloids"