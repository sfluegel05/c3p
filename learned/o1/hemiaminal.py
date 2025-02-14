"""
Classifies: CHEBI:73080 hemiaminal
"""
from rdkit import Chem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule is a hemiaminal based on its SMILES string.
    A hemiaminal is an organic compound where a single sp3-hybridized carbon atom is attached to both
    an amino group and a hydroxy group, not part of a carbonyl group. Hemiaminals can be cyclic or acyclic
    and may involve secondary or tertiary amines.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a hemiaminal, False otherwise.
        str: Reason for classification.
    """

    # Parse the SMILES string to create an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define SMARTS pattern for a hemiaminal functional group
    # [C;!$(C=O);X4] - sp3 carbon not double-bonded to oxygen (non-carbonyl)
    # (-[O]) - single bond to oxygen
    # (-[N]) - single bond to nitrogen
    hemiaminal_pattern = Chem.MolFromSmarts('[C;!$(C=O);X4](-[O])(-[N])')

    if hemiaminal_pattern is None:
        return False, "Error in SMARTS pattern."

    # Search for the hemiaminal pattern in the molecule
    matches = mol.GetSubstructMatches(hemiaminal_pattern)
    if matches:
        # For debugging or detailed information, you can uncomment the following lines
        # matched_atoms = [mol.GetAtomWithIdx(match[0]) for match in matches]
        # atom_indices = [atom.GetIdx() for atom in matched_atoms]
        return True, "Contains a carbon atom attached to both hydroxyl and amino groups (hemiaminal)."
    else:
        return False, "Does not contain a hemiaminal functional group."