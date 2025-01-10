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

    # Uracil base
    uracil_smarts = "O=C1NC=CC(=O)N1"
    uracil_pattern = Chem.MolFromSmarts(uracil_smarts)
    if uracil_pattern is None:
        return False, "Unable to create uracil pattern"

    # Ribose sugar attached to uracil
    ribose_smarts = "O[C@H]1[C@@H](O)[C@H](O)[C@@H](CO)O1"
    ribose_pattern = Chem.MolFromSmarts(ribose_smarts)
    if ribose_pattern is None:
        return False, "Unable to create ribose pattern"

    # Uridine (uracil + ribose)
    uridine_smarts = "O=C1NC=CC(=O)N1C2OC(CO)C(O)C2O"
    uridine_pattern = Chem.MolFromSmarts(uridine_smarts)
    if uridine_pattern is None:
        return False, "Unable to create uridine pattern"

    # Diphosphate group attached to ribose
    diphosphate_smarts = "OP(=O)(O)OP(=O)(O)O"
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)
    if diphosphate_pattern is None:
        return False, "Unable to create diphosphate pattern"

    # Sugar moiety attached via diphosphate linkage
    sugar_diphosphate_smarts = "OP(=O)([O-])OP(=O)([O-])OC1CO[C@@H](O)[C@H](O)[C@@H]1O"
    sugar_diphosphate_pattern = Chem.MolFromSmarts(sugar_diphosphate_smarts)
    if sugar_diphosphate_pattern is None:
        return False, "Unable to create sugar diphosphate pattern"

    # Check for uridine moiety (uracil + ribose)
    if not mol.HasSubstructMatch(uridine_pattern):
        return False, "Uridine moiety not found"

    # Check for diphosphate group attached to ribose
    diphosphate_matches = mol.GetSubstructMatches(diphosphate_pattern)
    if not diphosphate_matches:
        return False, "Diphosphate group not found attached to ribose"

    # Check for sugar moiety attached via diphosphate linkage
    # For simplicity, we check for any sugar ring attached to diphosphate
    sugar_ring_smarts = "OC1COC(O)C(O)C1O"
    sugar_ring_pattern = Chem.MolFromSmarts(sugar_ring_smarts)
    if sugar_ring_pattern is None:
        return False, "Unable to create sugar ring pattern"

    # Find atoms involved in diphosphate group
    diphosphate_atoms = [atom_idx for match in diphosphate_matches for atom_idx in match]

    # Check for sugar attached to diphosphate
    attached_sugars = False
    for phosphorous_idx in [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]:
        phosphorous_atom = mol.GetAtomWithIdx(phosphorous_idx)
        for neighbor in phosphorous_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for nbr in neighbor.GetNeighbors():
                    if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in diphosphate_atoms:
                        # Potential sugar attachment
                        if mol.HasSubstructMatch(sugar_ring_pattern, atomsToMatch=[nbr.GetIdx()]):
                            attached_sugars = True
                            break

    if not attached_sugars:
        return False, "No sugar moiety attached via diphosphate linkage"

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