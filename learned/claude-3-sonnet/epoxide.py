"""
Classifies: CHEBI:32955 epoxide
"""
"""
Classifies: CHEBI:35762 epoxide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_epoxide(smiles: str):
    """
    Determines if a molecule contains an epoxide group based on its SMILES string.
    An epoxide is a cyclic ether where the oxygen atom forms part of a 3-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an epoxide group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Try to sanitize the molecule
    try:
        Chem.SanitizeMol(mol)
    except:
        return False, "Could not sanitize molecule"

    # SMARTS pattern for epoxide - a three-membered ring containing an oxygen
    # [OX2r3] - oxygen with 2 connections in a 3-membered ring
    # 1 - first connection point
    # [#6r3] - any carbon in a 3-membered ring
    # [#6r3] - second carbon in a 3-membered ring
    # The whole pattern must form a 3-membered ring
    epoxide_pattern = Chem.MolFromSmarts('[OX2r3]1[#6r3][#6r3]1')
    
    if epoxide_pattern is None:
        return False, "Error in SMARTS pattern"
    
    # Look for matches with stereochemistry consideration
    matches = mol.GetSubstructMatches(epoxide_pattern, useChirality=True)
    
    if not matches:
        return False, "No epoxide group found"
        
    # Verify each match is actually a 3-membered ring containing an oxygen
    ring_info = mol.GetRingInfo()
    for match in matches:
        # Get the ring size for these atoms
        rings = ring_info.AtomRings()
        for ring in rings:
            if all(idx in ring for idx in match):
                if len(ring) == 3:  # Confirm it's a 3-membered ring
                    # Verify one atom is oxygen and others are carbons
                    atoms = [mol.GetAtomWithIdx(idx) for idx in match]
                    atom_types = [atom.GetAtomicNum() for atom in atoms]
                    if atom_types.count(8) == 1:  # One oxygen
                        return True, "Contains epoxide group (3-membered ring with oxygen)"
    
    return False, "No valid epoxide group found"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35762',
        'name': 'epoxide',
        'definition': 'Any cyclic ether in which the oxygen atom forms part of a 3-membered ring.',
        'parents': ['CHEBI:33641']
    }
}