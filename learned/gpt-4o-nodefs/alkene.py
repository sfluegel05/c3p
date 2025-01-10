"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene contains at least one carbon-carbon double bond typically in an open-chain structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkene, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define C=C double bond pattern
    alkene_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(alkene_pattern):
        return False, "No C=C double bond found"

    # Check for conjugated systems and ensure it is not polyene
    c_double_bonds = [bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE and 
                      bond.GetBeginAtom().GetAtomicNum() == 6 and bond.GetEndAtom().GetAtomicNum() == 6]
    if len(c_double_bonds) > 1:
        return False, "Contains multiple C=C bonds, possibly a polyene"

    # Get ring info
    ring_info = mol.GetRingInfo()
    small_ring = any(len(ring) < 8 for ring in ring_info.AtomRings())

    # Allow small rings only if they don't disrupt basic alkene structure
    if small_ring:
        # Pattern match to ignore non-disturbing ring structures like cyclohexene
        allowed_small_ring = Chem.MolFromSmarts("[C&R1]=[C&R1]")  # Simple ring alkene pattern
        if not mol.HasSubstructMatch(allowed_small_ring):
            return False, "Contains a small ring disrupting alkene structure"

    return True, "Contains at least one C=C double bond in a suitable mono-unsaturated structure"