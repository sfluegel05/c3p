"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:XX pyrroline
"""
from rdkit import Chem
from rdkit.Chem import rdChemReactions

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is any organic heteromonocyclic compound with a structure based on dihydropyrrole,
    which is a 5-membered ring containing one nitrogen atom and one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for pyrroline (dihydropyrrole ring)
    # 5-membered ring with one nitrogen atom and one double bond
    pyrroline_pattern = Chem.MolFromSmarts("[#6]-1-[#6]-[#6]-[#7]-[#6]=1")
    if pyrroline_pattern is None:
        return False, "Invalid SMARTS pattern for pyrroline"

    # Search for the pyrroline substructure in the molecule
    matches = mol.GetSubstructMatches(pyrroline_pattern)
    if matches:
        return True, "Contains pyrroline ring (5-membered ring with one nitrogen atom and one double bond)"
    else:
        return False, "Does not contain pyrroline ring"

    # Alternative method: Check all 5-membered rings
    """
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    for ring in atom_rings:
        if len(ring) == 5:
            num_nitrogen = 0
            num_double_bonds = 0
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    num_nitrogen += 1
            bonds = []
            for i in range(len(ring)):
                bond = mol.GetBondBetweenAtoms(ring[i], ring[(i+1)%len(ring)])
                if bond is not None:
                    bonds.append(bond)
            for bond in bonds:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    num_double_bonds += 1
            if num_nitrogen == 1 and num_double_bonds == 1:
                return True, "Contains pyrroline ring (5-membered ring with one nitrogen atom and one double bond)"
    return False, "Does not contain pyrroline ring"
    """
    
__metadata__ = {   'chemical_class': {   'id': 'CHEBI:XX',
                              'name': 'pyrroline',
                              'definition': 'Any organic heteromonocyclic compound with a structure based on a dihydropyrrole.'},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 0,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}