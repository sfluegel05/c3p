"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    Sesquiterpenoids are characterized by a 15-carbon skeleton derived from three isoprene units,
    and can have modifications such as removal of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 0. Basic carbon count
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count > 30: # reasonable upper bound
       return False, f"Too many carbons: {c_count}"

    # 1. Check for common sesquiterpenoid skeletons using SMARTS patterns.
    # These patterns cover several common backbones, including bicyclic and monocyclic
    core_patterns = [
        Chem.MolFromSmarts("[C;R0]([C;R0])([C;R0])([C;R0])[C;R0][C;R0]([C;R0])([C;R0])([C;R0])[C;R0][C;R0]([C;R0])([C;R0])[C;R0]"), # acyclic, 15 carbons
        Chem.MolFromSmarts("[C;R]([C;R])[C;R]([C;R])[C;R]1[C;R]([C;R])[C;R]([C;R])C([C;R])[C;R]1"), # simple monocyclic
        Chem.MolFromSmarts("[C;R]1[C;R]([C;R])[C;R]2[C;R]([C;R])[C;R]1[C;R]([C;R])[C;R]2"),  # bicyclic fused 5/5
        Chem.MolFromSmarts("[C;R]1[C;R]([C;R])[C;R]2[C;R]([C;R])[C;R]1[C;R]([C;R])[C;R]C2"), # bicyclic fused 5/6
        Chem.MolFromSmarts("[C;R]1[C;R]2[C;R]([C;R])[C;R]([C;R])[C;R]1[C;R]([C;R])[C;R]2"), # bicyclic fused 6/6
        Chem.MolFromSmarts("[C;R]1[C;R]2[C;R]3[C;R]1[C;R]2[C;R]3"), # tricyclic
    ]
    
    found_core = False
    for pattern in core_patterns:
        if mol.HasSubstructMatch(pattern):
            found_core = True
            break
    if not found_core:
            return False, "No common sesquiterpenoid skeleton found"

    # 2. Check ring system details (1-3 rings, mostly 5 and 6 members)
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if not (1 <= num_rings <= 3):
        return False, f"Incorrect number of rings: {num_rings} (should be between 1 and 3)."

    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not all(size in [5, 6] for size in ring_sizes):
      return False, f"Ring sizes not typical for sesquiterpenoid.  Ring sizes {ring_sizes}."

    # 3. Check for "isoprenoid" like chains (more flexible than previous attempt)
    #   Looking for a carbon chain with branches (allow for carbons to be CH or CH2 or C)
    chain_pattern = Chem.MolFromSmarts("[CX4]([CX4])~[CX4]~[CX4]([CX4])~[CX4]~[CX4]")
    if not mol.HasSubstructMatch(chain_pattern):
      return False, "Not branched as expected for isoprene units"

    return True, "Contains a sesquiterpenoid core, with typical ring system and branched chains."