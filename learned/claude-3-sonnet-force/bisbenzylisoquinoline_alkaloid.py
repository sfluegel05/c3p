"""
Classifies: CHEBI:133004 bisbenzylisoquinoline alkaloid
"""
"""
Classifies: bisbenzylisoquinoline alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_bisbenzylisoquinoline_alkaloid(smiles: str):
    """
    Determines if a molecule is a bisbenzylisoquinoline alkaloid based on its SMILES string.
    These compounds contain two benzylisoquinoline units linked by ether bridges.

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

    # Basic molecular properties check
    rings = rdMolDescriptors.CalcNumRings(mol)
    if rings < 6:
        return False, f"Too few rings ({rings}), need at least 6"

    # Count N and O atoms
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if n_count < 2:
        return False, f"Need at least 2 nitrogen atoms, found {n_count}"
    if o_count < 1:
        return False, f"Need oxygen atoms for ether bridges, found {o_count}"

    # More flexible patterns for isoquinoline/tetrahydroisoquinoline units
    # Matches both saturated and unsaturated variants with possible N-methylation
    isoquinoline_patterns = [
        # Basic tetrahydroisoquinoline core with optional N-methylation
        "[#6]1[#6][#6][N,N-C][#6][#6]2[#6](=[#6][#6]=[#6][#6]=2)[#6]1",
        # Alternative connection pattern
        "[#6]1[#6][#6][N,N-C][#6][#6]2[#6][#6]=[#6][#6]=[#6]12",
        # Pattern allowing for various substitutions
        "[#6]1[#6][#6][N,N-C][#6][#6]2[#6][#6][#6][#6][#6]12"
    ]

    # Convert patterns to RDKit molecules
    pattern_mols = [Chem.MolFromSmarts(p) for p in isoquinoline_patterns]
    
    # Count total matches across all patterns
    total_iso_units = 0
    for pattern in pattern_mols:
        if pattern is not None:
            matches = len(mol.GetSubstructMatches(pattern))
            total_iso_units += matches

    if total_iso_units < 2:
        return False, f"Need at least 2 isoquinoline-like units, found {total_iso_units}"

    # Look for benzyl groups attached to the isoquinoline cores
    benzyl_pattern = Chem.MolFromSmarts("[CH2]c1ccccc1")
    if benzyl_pattern:
        benzyl_count = len(mol.GetSubstructMatches(benzyl_pattern))
        if benzyl_count < 1:
            return False, "Missing benzyl groups"

    # Look for ether bridges between aromatic rings
    ether_bridge_patterns = [
        "c-O-c",  # Direct ether bridge
        "c-O-C-c", # Extended ether bridge
        "c-O-C(-*)c" # Branched ether bridge
    ]
    
    total_ether_bridges = 0
    for pattern in ether_bridge_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            total_ether_bridges += len(mol.GetSubstructMatches(pat))

    if total_ether_bridges < 1:
        return False, "No suitable ether bridges found"

    # Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for bisbenzylisoquinoline alkaloid"

    # Check for typical substituents
    substituent_patterns = {
        "methoxy": "cOC",
        "hydroxy": "cO[H]",
        "methylenedioxy": "OCO"
    }
    
    substituent_counts = {}
    for name, pattern in substituent_patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            substituent_counts[name] = len(mol.GetSubstructMatches(pat))
    
    total_substituents = sum(substituent_counts.values())
    if total_substituents < 1:
        return False, "Missing typical oxygen-containing substituents"

    reason = (f"Contains {total_iso_units} isoquinoline-like units, {total_ether_bridges} ether bridges, "
             f"and appropriate substituents. MW: {mol_wt:.1f}")
    
    return True, reason