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

    # SMARTS patterns for 1,4-naphthoquinone and 1,2-naphthoquinone cores
    # Note: Using more flexible patterns that account for both aromatic and non-aromatic bonds
    naphthoquinone_14 = Chem.MolFromSmarts(
        "[#6]1[#6](=[O,N])[#6][#6][#6]2[#6][#6][#6][#6][#6]2[#6]1=[O,N]"  # 1,4-naphthoquinone
    )
    naphthoquinone_12 = Chem.MolFromSmarts(
        "[#6]1[#6](=[O,N])[#6](=[O,N])[#6][#6]2[#6][#6][#6][#6][#6]2[#6]1"  # 1,2-naphthoquinone
    )

    # Check for either 1,4 or 1,2 naphthoquinone core
    has_14_core = mol.HasSubstructMatch(naphthoquinone_14)
    has_12_core = mol.HasSubstructMatch(naphthoquinone_12)
    
    if not (has_14_core or has_12_core):
        return False, "No naphthoquinone core found"

    # Get the core atoms
    core_matches = []
    if has_14_core:
        core_matches.extend(mol.GetSubstructMatches(naphthoquinone_14))
    if has_12_core:
        core_matches.extend(mol.GetSubstructMatches(naphthoquinone_12))
    
    if not core_matches:
        return False, "Could not identify naphthoquinone core atoms"
    
    # Create set of core atoms from first match
    core_atoms = set(core_matches[0])
    
    # Find hydroxy groups
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    
    if not hydroxy_matches:
        return False, "No hydroxy groups found"

    # Check if any hydroxy group is attached to the naphthoquinone core
    for hydroxy_match in hydroxy_matches:
        oh_oxygen = hydroxy_match[0]
        # Get the atom the OH is attached to
        oh_atom = mol.GetAtomWithIdx(oh_oxygen)
        for neighbor in oh_atom.GetNeighbors():
            # Check if the carbon the OH is attached to is part of the core
            # or is directly attached to the core
            if neighbor.GetIdx() in core_atoms:
                if has_14_core:
                    return True, "Found 1,4-hydroxynaphthoquinone structure"
                else:
                    return True, "Found 1,2-hydroxynaphthoquinone structure"
            # Check second-degree connections
            for next_neighbor in neighbor.GetNeighbors():
                if next_neighbor.GetIdx() in core_atoms:
                    if has_14_core:
                        return True, "Found 1,4-hydroxynaphthoquinone structure"
                    else:
                        return True, "Found 1,2-hydroxynaphthoquinone structure"

    return False, "Hydroxy groups present but not attached to naphthoquinone core"