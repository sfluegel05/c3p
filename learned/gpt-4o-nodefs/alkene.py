"""
Classifies: CHEBI:32878 alkene
"""
from rdkit import Chem

def is_alkene(smiles: str):
    """
    Determines if a molecule is an alkene based on its SMILES string.
    An alkene contains at least one carbon-carbon double bond in an open-chain structure,
    accounting for common substituents and functional groups that are permissible.

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

    # Find carbon-carbon double bond pattern (C=C)
    alkene_pattern = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(alkene_pattern):
        return False, "No C=C double bond found"

    # Identify and count SMARTS matches not to specifically exclude epoxides and allow simple alkene chains
    complex_subs_pattern = Chem.MolFromSmarts("[#7,#8,#9,#15,#16,#17,#35,#53]") # Heteroatoms: N, O, F, P, S, Cl, Br, I
    excluded_functionals = mol.HasSubstructMatch(complex_subs_pattern)

    # Get ring info
    ring_info = mol.GetRingInfo()
    small_ring = any(len(ring) < 8 for ring in ring_info.AtomRings())

    # Ensure any small ring found is an epoxide, otherwise reject
    if small_ring:
        epoxide_pattern = Chem.MolFromSmarts("[C;R2]1[O;R][C;R1]1")
        is_epoxide = mol.HasSubstructMatch(epoxide_pattern)
        if not is_epoxide:
            return False, "Contains a non-epoxide small ring"

    # Exclude due to complex substituents if true
    if excluded_functionals:
        return False, "Contains unacceptable heteroatoms or groups"

    return True, "Contains at least one C=C double bond in suitable alkene structure"