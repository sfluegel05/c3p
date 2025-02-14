"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: Isoflavones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_isoflavones(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.
    An isoflavone has a 3-aryl-1-benzopyran-4-one skeleton (3-aryl-4H-chromen-4-one)
    and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the isoflavone core
    # Chromen-4-one core: c1cc2ccccc2oc1=O
    # Aryl group at position 3
    isoflavone_pattern = Chem.MolFromSmarts('O=C1Oc2ccccc2cc1-c1ccccc1')
    if isoflavone_pattern is None:
        return False, "Error in SMARTS pattern"

    # Check if the molecule matches the isoflavone core pattern
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Does not contain the 3-aryl-1-benzopyran-4-one skeleton"

    # Verify the presence of the chromen-4-one core
    chromen4one_pattern = Chem.MolFromSmarts('O=C1Oc2ccccc2cc1')
    if not mol.HasSubstructMatch(chromen4one_pattern):
        return False, "Does not contain the chromen-4-one core"

    # Verify the aryl group attached at position 3
    aryl_group_pattern = Chem.MolFromSmarts('c1ccccc1')
    matches = mol.GetSubstructMatches(chromen4one_pattern)
    if matches:
        for match in matches:
            atom_indices = set(match)
            # Check for aryl substituent at position 3 of the chromen-4-one core
            position_3_atom = mol.GetAtomWithIdx(match[5])  # Atom at position 3
            neighbors = position_3_atom.GetNeighbors()
            for neighbor in neighbors:
                if neighbor.GetIdx() not in atom_indices:
                    # Check if the substituent is an aryl group
                    neighbor_submol = Chem.PathToSubmol(mol, [position_3_atom.GetIdx(), neighbor.GetIdx()])
                    if neighbor_submol.HasSubstructMatch(aryl_group_pattern):
                        return True, "Contains the 3-aryl-1-benzopyran-4-one skeleton with an aryl group at position 3"
        return False, "No aryl group attached at position 3 of the chromen-4-one core"
    else:
        return False, "Chromen-4-one core not found"

    return False, "Molecule does not match isoflavone criteria"