"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol must have at least one hydroxyl group and one nitro group attached to an aromatic system.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string to a molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"

    # Check for phenolic hydroxyl in an aromatic ring
    hydroxyl_pattern = Chem.MolFromSmarts("[$([OH])][c]")  # This pattern ensures the hydroxyl is connected to an aromatic carbon but more loosely
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No phenolic hydroxyl group found on aromatic system"

    # Check for a nitro group
    nitro_group_pattern = Chem.MolFromSmarts("[NX3](=O)[O-]")  # Standard nitro group pattern
    if not mol.HasSubstructMatch(nitro_group_pattern):
        return False, "No nitro group found, need at least one nitro group"

    # Ensure that hydroxyl and nitro groups are on the same benzene-type aromatic ring
    atom_indices = [atom.GetIdx() for atom in mol.GetAtoms()]
    hydroxyl_atoms = mol.GetSubstructMatches(hydroxyl_pattern)
    nitro_atoms = mol.GetSubstructMatches(nitro_group_pattern)
    
    for hydroxyl in hydroxyl_atoms:
        for nitro in nitro_atoms:
            if any(Chem.FindAtomEnvironmentOfRadiusN(mol, 2, hydroxyl[0], nitro[0])):
                return True, "Phenolic hydroxyl and nitro group on same aromatic system found"

    return False, "Phenolic hydroxyl and nitro group not on the same benzene ring"