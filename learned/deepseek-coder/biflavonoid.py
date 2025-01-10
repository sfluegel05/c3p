"""
Classifies: CHEBI:50128 biflavonoid
"""
"""
Classifies: CHEBI:77663 biflavonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid oligomer obtained by the oxidative coupling of at least two units of aryl-substituted benzopyran rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more flexible flavonoid pattern (benzopyran ring system with possible substitutions)
    flavonoid_pattern = Chem.MolFromSmarts("[c]1[c][c][c][c][c]1-[O][C]2[C]=[C][C]=[C][C]=2")
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    
    if len(flavonoid_matches) < 2:
        return False, f"Found {len(flavonoid_matches)} flavonoid units, need at least 2"

    # Check if the flavonoid units are connected by a single bond or atom
    connected = False
    for i in range(len(flavonoid_matches)):
        for j in range(i+1, len(flavonoid_matches)):
            # Check if there is a bond between any atom in the first unit and any atom in the second unit
            for atom1 in flavonoid_matches[i]:
                for atom2 in flavonoid_matches[j]:
                    if mol.GetBondBetweenAtoms(atom1, atom2) is not None:
                        connected = True
                        break
                if connected:
                    break
            if connected:
                break
        if connected:
            break

    if not connected:
        return False, "Flavonoid units are not connected by a single bond or atom"

    # Check if the molecule is an oligomer (i.e., it consists of two flavonoid units)
    # This is a heuristic check based on the number of flavonoid units
    if len(flavonoid_matches) != 2:
        return False, f"Found {len(flavonoid_matches)} flavonoid units, need exactly 2 for biflavonoid"

    return True, "Contains two flavonoid units connected by a single bond or atom"