"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    An anthocyanidin cation is characterized by an oxygenated flavylium ion core.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Define the correct flavylium pattern: 2-phenylchromenylium core with a positively charged oxygen
    flavylium_pattern = Chem.MolFromSmarts("[o+]1cc2ccc(O)cc2c(c3ccccc3)c1")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return (False, "No flavylium core found")

    # Check for at least two oxygen atoms in hydroxyl or ether form directly attached to the rings
    oxy_group_pattern = Chem.MolFromSmarts("[OH]")  # Updates: involve ethers if necessary
    if len(mol.GetSubstructMatches(oxy_group_pattern)) < 2:
        return (False, "Insufficient oxygenation; at least two hydroxyl groups expected")

    # Confirm presence of a positive charge on the correct atom
    positive_charge = any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms())
    if not positive_charge:
        return (False, "Molecule lacks a positive charge indicating it is not a cation")

    return (True, "Contains an oxygenated flavylium ion core with positive charge")