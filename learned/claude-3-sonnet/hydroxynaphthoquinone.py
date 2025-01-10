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

    # First check for naphthalene core (two fused 6-membered rings)
    naphthalene_pattern = Chem.MolFromSmarts("c1ccc2ccccc2c1")
    if not mol.HasSubstructMatch(naphthalene_pattern):
        return False, "No naphthalene core found"

    # Check for ketone groups (=O)
    ketone_pattern = Chem.MolFromSmarts("C(=O)")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if len(ketone_matches) < 2:
        return False, f"Found only {len(ketone_matches)} ketone groups, need at least 2"

    # Check for hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if not hydroxy_matches:
        return False, "No hydroxy groups found"

    # Get the naphthalene core atoms
    core_matches = mol.GetSubstructMatches(naphthalene_pattern)
    core_atoms = set(core_matches[0])
    
    # Find atoms that are part of the extended naphthoquinone system
    # (includes the core and directly attached atoms)
    extended_core = set()
    for atom_idx in core_atoms:
        extended_core.add(atom_idx)
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            extended_core.add(neighbor.GetIdx())

    # Check if ketones are part of the naphthoquinone system
    ketone_on_core = False
    for ketone_match in ketone_matches:
        if ketone_match[0] in extended_core:
            ketone_on_core = True
            break
    
    if not ketone_on_core:
        return False, "Ketone groups not part of naphthoquinone system"

    # Check if any hydroxy group is attached to the extended naphthoquinone core
    for hydroxy_match in hydroxy_matches:
        oh_oxygen = hydroxy_match[0]
        # Get the atom the OH is attached to
        neighbors = mol.GetAtomWithIdx(oh_oxygen).GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIdx() in extended_core:
                return True, "Contains naphthoquinone core with at least one hydroxy group attached"

    return False, "Hydroxy groups present but not attached to naphthoquinone core"