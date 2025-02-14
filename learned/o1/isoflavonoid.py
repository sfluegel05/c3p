"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
"""

from rdkit import Chem
from rdkit.Chem import rdqueries

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is any molecule containing the isoflavonoid core structure,
    which is a 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the benzopyran ring (1-benzopyran)
    benzopyran_smarts = 'c1ccc2occ(c3ccccc13)cc2'
    benzopyran_pattern = Chem.MolFromSmarts(benzopyran_smarts)
    if benzopyran_pattern is None:
        return False, "Invalid benzopyran SMARTS pattern"

    # Search for benzopyran core
    matches = mol.GetSubstructMatches(benzopyran_pattern)
    if not matches:
        return False, "Does not contain 1-benzopyran core structure"

    # For each benzopyran match, check for aryl substituent at position 3
    for match in matches:
        # Map the atom indices
        atom_indices = {atom.GetAtomMapNum(): atom.GetIdx() for atom in benzopyran_pattern.GetAtoms()}
        benzopyran_atoms = [match[i] for i in range(len(match))]
        
        # Position 3 in benzopyran SMARTS is atom with AtomMapNum 10
        position_3 = match[9]  # Adjusted index for zero-based indexing

        # Get the neighbors of position 3 atom excluding those in the benzopyran ring
        position_3_atom = mol.GetAtomWithIdx(position_3)
        neighbors = position_3_atom.GetNeighbors()
        for neighbor in neighbors:
            if neighbor.GetIdx() not in benzopyran_atoms:
                # Check if the substituent is an aryl group (aromatic ring)
                if neighbor.GetIsAromatic():
                    # Check if the substituent is an aromatic ring
                    ring_info = mol.GetRingInfo()
                    if ring_info.IsAtomInRingOfSize(neighbor.GetIdx(), 6):
                        return True, "Contains isoflavonoid core with aryl substituent at position 3"
    return False, "Does not have aryl substituent at position 3 of benzopyran core"