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

    # Look for tetrahydroisoquinoline cores (can be N-methylated or not)
    # Pattern matches both aromatic and saturated versions
    thiq_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#7][#6][#6]2[#6]([#6][#6]=1)[#6][#6][#6]=2")
    thiq_matches = mol.GetSubstructMatches(thiq_pattern)
    
    if len(thiq_matches) < 2:
        return False, "Must contain two tetrahydroisoquinoline cores"

    # Look for ether bridges (-O-) connecting aromatic rings
    ether_bridge_pattern = Chem.MolFromSmarts("c-[#8]-c")
    ether_matches = mol.GetSubstructMatches(ether_bridge_pattern)
    
    if len(ether_matches) < 1:
        return False, "Must contain ether bridges between aromatic rings"

    # Check for nitrogen-containing 6-membered rings
    n_ring_pattern = Chem.MolFromSmarts("[#7]1[#6][#6][#6][#6][#6]1")
    n_ring_matches = mol.GetSubstructMatches(n_ring_pattern)
    
    if len(n_ring_matches) < 2:
        return False, "Must contain two nitrogen-containing 6-membered rings"

    # Check for benzyl groups
    benzyl_pattern = Chem.MolFromSmarts("c1ccccc1-[CH2]")
    benzyl_matches = mol.GetSubstructMatches(benzyl_pattern)
    
    if len(benzyl_matches) < 2:
        return False, "Must contain two benzyl groups"

    # Check for common substituents (methoxy, hydroxy)
    methoxy_pattern = Chem.MolFromSmarts("c-[#8]-[CH3]")
    hydroxy_pattern = Chem.MolFromSmarts("c-[#8;H1]")
    substituents = len(mol.GetSubstructMatches(methoxy_pattern)) + len(mol.GetSubstructMatches(hydroxy_pattern))
    
    if substituents < 2:
        return False, "Insufficient number of typical substituents (methoxy/hydroxy)"

    # Count aromatic rings
    aromatic_rings = 0
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
            aromatic_rings += 1
            
    if aromatic_rings < 4:  # Should have at least 4 aromatic rings
        return False, "Insufficient number of aromatic rings"

    # Verify overall molecular properties
    if not (400 <= mol.GetExactMolWt() <= 800):  # Typical mass range for this class
        return False, "Molecular weight outside typical range"

    # Additional structural checks
    # Look for characteristic carbon skeleton connecting the two units
    connecting_pattern = Chem.MolFromSmarts("c1cccc(Oc2ccccc2)c1")
    if not mol.HasSubstructMatch(connecting_pattern):
        return False, "Missing characteristic biaryl ether linkage"

    return True, "Contains two benzylisoquinoline units connected by ether bridges with appropriate substitution patterns"