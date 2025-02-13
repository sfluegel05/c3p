"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: CHEBI:26421 quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for cyclic ketone patterns
    # This pattern matches both para and ortho quinones
    para_quinone = Chem.MolFromSmarts("[#6;R]=O.[#6;R]=O")
    ortho_quinone = Chem.MolFromSmarts("[#6;R](=O)-[#6;R]-[#6;R](=O)")
    
    has_para = mol.HasSubstructMatch(para_quinone)
    has_ortho = mol.HasSubstructMatch(ortho_quinone)
    
    if not (has_para or has_ortho):
        return False, "Must contain cyclic ketone groups in quinone arrangement"

    # Count ring systems containing ketones
    ketone_pattern = Chem.MolFromSmarts("[#6;R]=O")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    
    if len(ketone_matches) < 2:
        return False, "Must have at least two ketone groups"
    
    # Check if ketones are in conjugated system
    ring_info = mol.GetRingInfo()
    ring_systems = ring_info.AtomRings()
    
    # Get the atoms between ketones
    connecting_atoms = []
    for ring in ring_systems:
        ketones_in_ring = []
        for ketone in ketone_matches:
            if ketone[0] in ring:
                ketones_in_ring.append(ketone[0])
        if len(ketones_in_ring) >= 2:
            connecting_atoms.extend(ring)
    
    if not connecting_atoms:
        return False, "Ketone groups must be in same ring system"
    
    # Check for conjugation or aromaticity in connecting atoms
    has_conjugation = False
    for atom_idx in connecting_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetIsAromatic():
            has_conjugation = True
            break
        for bond in atom.GetBonds():
            if bond.GetIsConjugated():
                has_conjugation = True
                break
    
    if not has_conjugation:
        return False, "Must have conjugated system between ketones"
    
    # Additional check for polycyclic systems
    polycyclic = len(ring_systems) > 1
    
    # Construct reason string
    if has_para:
        structure_type = "para"
    else:
        structure_type = "ortho"
    
    system_type = "polycyclic" if polycyclic else "monocyclic"
    reason = f"Contains {structure_type}-quinone structure in {system_type} conjugated system"
    
    return True, reason