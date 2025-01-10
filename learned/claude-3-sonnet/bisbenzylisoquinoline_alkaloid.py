"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: CHEBI:bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a bisbenzylisoquinoline alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for minimum complexity
    if mol.GetNumAtoms() < 30:
        return False, "Molecule too small for bisbenzylisoquinoline alkaloid"

    # More flexible patterns for isoquinoline-like cores
    # Basic tetrahydroisoquinoline pattern with variations
    patterns = [
        # Basic tetrahydroisoquinoline core with flexible bonds
        "[#6]1~[#6]~[#6]~[#7]~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~1~[#6]2",
        # Pattern with N-methyl variation
        "[#6]1~[#6]~[#6]~[#7X3]($([CH3]))~[#6]~[#6]2~[#6]~[#6]~[#6]~[#6]~1~[#6]2",
        # Pattern allowing for different oxidation states
        "[#6]1~[#6]~[#6]~[#7]~[#6]~[#6]2~[#6]=,:[#6]~[#6]~[#6]~1~[#6]2"
    ]
    
    total_matches = 0
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            matches = len(mol.GetSubstructMatches(patt))
            total_matches += matches

    if total_matches < 2:
        return False, "Must contain at least two isoquinoline-like cores"

    # Check for ether bridges between aromatic rings (more flexible pattern)
    ether_patterns = [
        "c-[#8]-c",  # Basic ether bridge
        "c-[#8]-[#6]-c",  # Extended ether bridge
        "c1ccc2OCOc2c1"   # Methylenedioxy bridge
    ]
    
    ether_matches = 0
    for pattern in ether_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            ether_matches += len(mol.GetSubstructMatches(patt))

    if ether_matches < 1:
        return False, "Must contain ether bridges between units"

    # Check for typical substituents with more flexible patterns
    substituent_patterns = [
        "c-[#8]-[CH3]",  # Methoxy
        "c-[#8;H1]",     # Hydroxy
        "c1ccc2OCOc2c1"  # Methylenedioxy
    ]
    
    total_substituents = 0
    for pattern in substituent_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt:
            total_substituents += len(mol.GetSubstructMatches(patt))

    if total_substituents < 2:
        return False, "Insufficient number of characteristic substituents"

    # Count nitrogen atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, "Must contain at least two nitrogen atoms"

    # Count aromatic rings
    ring_info = mol.GetRingInfo()
    aromatic_rings = 0
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings += 1

    if aromatic_rings < 4:
        return False, "Insufficient number of aromatic rings"

    # Check molecular weight range
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if not (400 <= mol_wt <= 800):
        return False, "Molecular weight outside typical range"

    # Verify molecule is a single connected system
    if len(Chem.GetMolFrags(mol)) > 1:
        return False, "Structure appears to be fragmented"

    return True, "Contains two benzylisoquinoline units with appropriate bridging and substitution patterns"