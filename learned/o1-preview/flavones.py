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

    # Define the flavone core using a SMARTS pattern with aromaticity
    flavone_core_smarts = 'O=C1c2ccccc2Oc2ccccc12'  # Aromatic flavone core
    flavone_core = Chem.MolFromSmarts(flavone_core_smarts)
    if flavone_core is None:
        return False, "Error in flavone core SMARTS definition"

    # Find substructure matches of the flavone core
    core_matches = mol.GetSubstructMatches(flavone_core)
    if not core_matches:
        return False, "Flavone core structure not found"

    # Optionally, check for substitutions at position 2 (the aryl group)
    # Position 2 corresponds to atom index 3 in the flavone core
    for match in core_matches:
        position_2_idx = match[3]  # Atom index of position 2 in the molecule
        position_2_atom = mol.GetAtomWithIdx(position_2_idx)
        # Check if position 2 has a substituent that is an aryl group
        for bond in position_2_atom.GetBonds():
            neighbor = bond.GetOtherAtom(position_2_atom)
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in match:
                # Check if the substituent is an aryl group (aromatic ring)
                env = Chem.FindAtomEnvironmentOfRadiusN(mol, radius=2, atomId=neighbor_idx)
                amap = {}
                submol = Chem.PathToSubmol(mol, env, atomMap=amap)
                aryl_pattern = Chem.MolFromSmarts('a1aaaaa1')  # Aromatic ring of size 6
                if submol.HasSubstructMatch(aryl_pattern):
                    return True, "Contains flavone core with 2-aryl substitution"
        # If no aryl substituent found at position 2
        return False, "No aryl substitution at position 2 of flavone core"

    # If no matches meet the criteria
    return False, "Flavone core with 2-aryl substitution not found"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'flavones',
        'definition': 'A member of the class of flavonoid with a 2-aryl-1-benzopyran-4-one (2-arylchromen-4-one) skeleton and its substituted derivatives.',
        'parents': []
    }
}