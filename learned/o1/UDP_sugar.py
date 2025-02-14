"""
Classifies: CHEBI:17297 UDP-sugar
"""
"""
Classifies: UDP-sugar
"""
from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    A UDP-sugar is a pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to a sugar via an anomeric diphosphate linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a UDP-sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for the components of UDP-sugar

    # Uracil base (without stereochemistry)
    uracil_smarts = "O=C1NC=CC(=O)N1"
    uracil_pattern = Chem.MolFromSmarts(uracil_smarts)

    # Ribose sugar attached to uracil (without stereochemistry)
    ribose_smarts = "C1OC(CO)C(O)C1O"
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)

    # Uridine (uracil base attached to ribose)
    uridine_smarts = "O=C1NC=CC(=O)N1C2OC(CO)C(O)C2O"
    uridine_pattern = Chem.MolFromSmarts(uridine_smarts)

    # Diphosphate group
    diphosphate_smarts = "OP(=O)(O)OP(=O)(O)O"
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)

    # Sugar moiety (any monosaccharide without stereochemistry)
    sugar_smarts = "C1OC(O)C(O)C(O)C1O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)

    # Check for uridine moiety
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "Uridine moiety not found"

    # Get the uridine match atoms
    uridine_match = mol.GetSubstructMatch(uridine_pattern)
    uridine_atoms = set(uridine_match)

    # Check for diphosphate group
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_matches:
        return False, "Diphosphate group not found"

    # For each diphosphate match, check if it's connected to uridine
    diphosphate_connected = False
    for diphosphate_match in diphosphate_matches:
        diphosphate_atoms = set(diphosphate_match)
        # Check for connection between uridine and diphosphate
        for atom_idx in uridine_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in diphosphate_atoms:
                    diphosphate_connected = True
                    break
            if diphosphate_connected:
                break
        if diphosphate_connected:
            break

    if not diphosphate_connected:
        return False, "Diphosphate group not attached to uridine"

    # Check for sugar moiety attached via diphosphate
    sugar_connected = False
    for diphosphate_match in diphosphate_matches:
        diphosphate_atoms = set(diphosphate_match)
        # Find atoms connected to diphosphate but not part of uridine
        for atom_idx in diphosphate_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx not in diphosphate_atoms and nbr_idx not in uridine_atoms:
                    # Potential sugar attachment
                    if mol.HasSubstructMatch(sugar_pattern, useChirality=False, atomIdxs=[nbr_idx]):
                        sugar_connected = True
                        break
            if sugar_connected:
                break
        if sugar_connected:
            break

    if not sugar_connected:
        return False, "Sugar moiety not attached via diphosphate linkage"

    return True, "UDP-sugar identified with uridine diphosphate linked to sugar via anomeric diphosphate"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:46211',
        'name': 'UDP-sugar',
        'definition': 'A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to an unspecified sugar via an anomeric diphosphate linkage.',
        'parents': []
    },
    'config': {}
}