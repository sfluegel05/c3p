"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule contains a guaiacol moiety based on its SMILES string.
    A guaiacol is a phenol with a methoxy group in the ortho position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains guaiacol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More specific SMARTS pattern for guaiacol core:
    # [cH,c]:1 - aromatic carbon that's part of a ring
    # ([OH1]) - hydroxyl group
    # (:c([OX2][CH3])) - ortho carbon with methoxy group
    # The rest just ensures it's part of a benzene ring
    guaiacol_pattern = Chem.MolFromSmarts('[cH,c](:c([OH1]):c(:c([OX2][CH3])):c:c:c:1)')
    alt_pattern = Chem.MolFromSmarts('[cH,c](:c([OX2][CH3]):c(:c([OH1])):c:c:c:1)')

    if not (mol.HasSubstructMatch(guaiacol_pattern) or mol.HasSubstructMatch(alt_pattern)):
        return False, "Does not contain guaiacol pattern"

    # Get matches
    matches = mol.GetSubstructMatches(guaiacol_pattern)
    alt_matches = mol.GetSubstructMatches(alt_pattern)
    all_matches = list(matches) + list(alt_matches)

    for match in all_matches:
        # Get the atoms involved in the match
        ring_atoms = set(match)
        
        # Verify these atoms are part of the same ring
        ri = mol.GetRingInfo()
        for ring in ri.AtomRings():
            if all(idx in ring for idx in ring_atoms):
                # Verify the ring is aromatic
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    # Verify methoxy is truly -OCH3
                    for atom_idx in match:
                        atom = mol.GetAtomWithIdx(atom_idx)
                        if atom.GetAtomicNum() == 8:  # Oxygen
                            for neighbor in atom.GetNeighbors():
                                if neighbor.GetAtomicNum() == 6:  # Carbon
                                    if neighbor.GetDegree() == 4:  # sp3 carbon
                                        return True, "Contains phenol with ortho methoxy group (guaiacol moiety)"

    return False, "Contains required elements but not in correct configuration"