"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
"""
Classifies: hydroxynaphthoquinone
Definition: Any naphthoquinone in which the naphthaoquinone moiety is substituted by at least one hydroxy group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_hydroxynaphthoquinone, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for naphthoquinone core - two fused 6-membered rings where some carbons can be non-aromatic
    # due to ketone groups
    naphthoquinone_pattern = Chem.MolFromSmarts("[#6]1[#6][#6][#6]2[#6][#6][#6][#6][#6][#6]2[#6]1")
    if not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No bicyclic core found"

    # Find naphthoquinone system with two ketones
    # This pattern looks for a bicyclic system with two ketone groups
    naphthoquinone_with_ketones = Chem.MolFromSmarts(
        "[#6]1[#6][#6][#6]2[#6][#6][#6][#6][#6][#6]2[#6]1(=[O,N])" # Bicyclic core with one ketone
        ".[#6]1[#6][#6][#6]2[#6][#6][#6][#6][#6][#6]2[#6]1(=[O,N])" # Same core with second ketone
    )
    if not mol.HasSubstructMatch(naphthoquinone_with_ketones):
        return False, "No naphthoquinone system found (missing ketone groups)"

    # Get the naphthoquinone core atoms
    core_matches = mol.GetSubstructMatches(naphthoquinone_pattern)
    if not core_matches:
        return False, "Could not identify naphthoquinone core atoms"
    
    core_atoms = set(core_matches[0])
    
    # Find atoms that are part of the extended naphthoquinone system
    extended_core = set()
    for atom_idx in core_atoms:
        extended_core.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            extended_core.add(neighbor.GetIdx())

    # Check for hydroxy groups attached to the naphthoquinone core
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    if not hydroxy_matches:
        return False, "No hydroxy groups found"

    # Check if any hydroxy group is attached to the extended naphthoquinone core
    for hydroxy_match in hydroxy_matches:
        oh_oxygen = hydroxy_match[0]
        # Get the atom the OH is attached to
        neighbors = mol.GetAtomWithIdx(oh_oxygen).GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIdx() in core_atoms:  # Check core atoms only, not extended core
                return True, "Found hydroxynaphthoquinone structure"

    return False, "Hydroxy groups present but not attached to naphthoquinone core"