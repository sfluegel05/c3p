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
    These compounds consist of two benzylisoquinoline units linked by ether bridges.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule matches criteria, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Adjusted isoquinoline pattern with possible substitutions
    isoquinoline_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6]2[#6](=[#7]-[#6])[#6][#6][#6][#6]12")
    if not isoquinoline_pattern:
        return None, "Invalid isoquinoline pattern"
    
    # Find at least two isoquinoline units
    iso_matches = mol.GetSubstructMatches(isoquinoline_pattern)
    if len(iso_matches) < 2:
        return False, f"Found {len(iso_matches)} isoquinoline units, need at least 2"
    
    # Check for ether bridges (O connecting two aromatic systems)
    ether_bridge_pattern = Chem.MolFromSmarts("[#6]~[#8]~[#6]")
    ether_matches = mol.GetSubstructMatches(ether_bridge_pattern)
    if len(ether_matches) < 1:
        return False, "No ether bridges found"
    
    # Verify tertiary amines (common in these alkaloids)
    tertiary_amine_pattern = Chem.MolFromSmarts("[#7X3]([#6])([#6])[#6]")
    amine_matches = mol.GetSubstructMatches(tertiary_amine_pattern)
    if len(amine_matches) < 2:
        return False, f"Found {len(amine_matches)} tertiary amines, need at least 2"
    
    # Check molecular complexity (minimum 40 atoms typical)
    atom_count = mol.GetNumAtoms()
    if atom_count < 40:
        return False, f"Only {atom_count} atoms, insufficient complexity"
    
    # Verify multiple aromatic rings (at least 4 benzene/isoquinoline rings)
    aromatic_rings = [ring for ring in mol.GetRingInfo().AtomRings() 
                     if all(mol.GetAtomWithIdx(a).GetIsAromatic() for a in ring)]
    if len(aromatic_rings) < 4:
        return False, f"Only {len(aromatic_rings)} aromatic rings found"
    
    return True, "Contains two isoquinoline units connected by ether bridges with tertiary amines"