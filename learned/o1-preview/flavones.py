"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: flavones
"""
from rdkit import Chem

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone based on its SMILES string.
    A flavone has a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton
    and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the chromen-4-one core (flavone core)
    flavone_core_smiles = 'O=C1C=CC2=CC=CC=C2O1'  # SMILES for chromen-4-one
    flavone_core = Chem.MolFromSmiles(flavone_core_smiles)
    if flavone_core is None:
        return False, "Error in flavone core definition"

    # Find substructure matches of the flavone core
    core_matches = mol.GetSubstructMatches(flavone_core)
    if not core_matches:
        return False, "Flavone core structure not found"

    # Define an aromatic ring pattern (aryl group)
    aryl_pattern = Chem.MolFromSmarts('a1aaaaa1')  # Aromatic ring of size 6

    # For each match, check for aryl substitution at position 2
    # Position 2 corresponds to atom index 2 in flavone_core
    for match in core_matches:
        # Get the atom index of position 2 in the molecule
        position_2_idx = match[2]
        position_2_atom = mol.GetAtomWithIdx(position_2_idx)
        # Check neighbors of position 2 atom
        for neighbor in position_2_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            # If the neighbor atom is not part of the core match, it's a substituent
            if neighbor_idx not in match:
                # Check if the substituent is an aryl group
                # Get the fragment connected to the neighbor atom
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius=3, atomId=neighbor_idx)
                amap = {}
                submol = Chem.PathToSubmol(mol, env, atomMap=amap)
                if submol.HasSubstructMatch(aryl_pattern):
                    return True, "Contains flavone core with 2-aryl substitution"
    # If no aryl group found at position 2
    return False, "No aryl substitution at position 2 of flavone core"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'flavones',
        'definition': 'A member of the class of flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton and its substituted derivatives.',
        'parents': []
    }
}