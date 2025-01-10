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
    if mol.GetNumAtoms() < 30:  # These are large molecules
        return False, "Molecule too small for bisbenzylisoquinoline alkaloid"

    # Look for isoquinoline-like cores (more flexible patterns)
    # Pattern 1: Basic tetrahydroisoquinoline core
    thiq_pattern1 = Chem.MolFromSmarts("[#6]1[#6][#6][#7][#6][#6]2[#6][#6][#6][#6]=1[#6]2")
    # Pattern 2: Alternative with different oxidation states
    thiq_pattern2 = Chem.MolFromSmarts("[#6]1[#6][#6][#7][#6][#6]2[#6][#6]=,:[#6][#6]=1[#6]2")
    # Pattern 3: N-methylated version
    thiq_pattern3 = Chem.MolFromSmarts("[#6]1[#6][#6][N;!H0,$(NC)][#6][#6]2[#6][#6][#6][#6]=1[#6]2")
    
    thiq_matches = (len(mol.GetSubstructMatches(thiq_pattern1)) + 
                   len(mol.GetSubstructMatches(thiq_pattern2)) +
                   len(mol.GetSubstructMatches(thiq_pattern3)))
    
    if thiq_matches < 2:
        return False, "Must contain two isoquinoline-like cores"

    # Look for benzyl groups attached to nitrogen (more flexible pattern)
    benzyl_n_pattern = Chem.MolFromSmarts("[#7]C[c]1[c][c][c][c][c]1")
    if len(mol.GetSubstructMatches(benzyl_n_pattern)) < 1:
        return False, "Must contain benzyl groups attached to nitrogen"

    # Look for ether bridges between aromatic rings (including various positions)
    ether_bridge_pattern = Chem.MolFromSmarts("c-[#8]-c")
    ether_matches = mol.GetSubstructMatches(ether_bridge_pattern)
    if len(ether_matches) < 1:
        return False, "Must contain ether bridges between aromatic rings"

    # Check for typical substituents (methoxy, hydroxy) on aromatic rings
    methoxy_pattern = Chem.MolFromSmarts("c-[#8]-[CH3]")
    hydroxy_pattern = Chem.MolFromSmarts("c-[#8;H1]")
    substituents = len(mol.GetSubstructMatches(methoxy_pattern)) + len(mol.GetSubstructMatches(hydroxy_pattern))
    
    if substituents < 2:
        return False, "Insufficient number of typical substituents (methoxy/hydroxy)"

    # Count nitrogen atoms (should have at least 2 for the two isoquinoline units)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 2:
        return False, "Must contain at least two nitrogen atoms"

    # Count aromatic rings
    aromatic_rings = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings += 1
            
    if aromatic_rings < 4:  # Should have at least 4 aromatic rings
        return False, "Insufficient number of aromatic rings"

    # Verify molecular weight is in typical range
    mol_wt = Chem.Descriptors.ExactMolWt(mol)
    if not (400 <= mol_wt <= 800):
        return False, "Molecular weight outside typical range"

    # Additional check for connected system
    # The molecule should be one connected system, not fragments
    if len(Chem.GetMolFrags(mol)) > 1:
        return False, "Structure appears to be fragmented"

    return True, "Contains two benzylisoquinoline units with appropriate bridging and substitution patterns"