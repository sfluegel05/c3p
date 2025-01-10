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

    # Multiple SMARTS patterns to catch different epoxide representations
    patterns = [
        # Basic 3-membered ring with oxygen
        "[O;R]1[C;R][C;R]1",
        # Alternative representation
        "C1OC1",
        # Catch spiro cases
        "[O;R]12[C]1[C]2"
    ]
    
    for pattern in patterns:
        epoxide_pattern = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(epoxide_pattern):
            matches = mol.GetSubstructMatches(epoxide_pattern)
            
            for match in matches:
                # Get the atoms in the potential epoxide
                atoms = [mol.GetAtomWithIdx(idx) for idx in match]
                
                # Verify it's a 3-membered ring
                ring_info = mol.GetRingInfo()
                
                # Check if these atoms form a 3-membered ring
                for ring in ring_info.AtomRings():
                    if all(atom.GetIdx() in ring for atom in atoms) and len(ring) == 3:
                        # Verify one atom is oxygen and others are carbons
                        atom_types = [atom.GetAtomicNum() for atom in atoms]
                        if 8 in atom_types:  # 8 is atomic number for oxygen
                            return True, "Contains epoxide group (3-membered ring with oxygen)"
    
    return False, "No epoxide group found (3-membered ring with oxygen)"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35762',
        'name': 'epoxide',
        'definition': 'Any cyclic ether in which the oxygen atom forms part of a 3-membered ring.',
        'parents': ['CHEBI:33641']
    }
}