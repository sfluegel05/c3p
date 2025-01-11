"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: CHEBI:35218 anthocyanidin cation
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidin cations are oxygenated derivatives of flavylium (2-phenylchromenylium) with a positive charge.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a positively charged oxygen in the chromenylium ring
    # The pattern matches the flavylium cation structure: [O+]=C1C=C(C2=CC=CC=C2)C=C1
    flavylium_pattern = Chem.MolFromSmarts("[O+]=c1ccc2ccccc2c1")
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium cation structure found"

    # Check for hydroxyl groups attached to the aromatic rings
    # The pattern matches at least one hydroxyl group attached to the aromatic rings
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "No hydroxyl groups found on aromatic rings"

    # Check for the presence of at least one oxygenated aromatic ring
    # The pattern matches an aromatic ring with at least one oxygen substituent
    oxygenated_aromatic_pattern = Chem.MolFromSmarts("c1ccccc1[OH]")
    if not mol.HasSubstructMatch(oxygenated_aromatic_pattern):
        return False, "No oxygenated aromatic ring found"

    # Check for the presence of at least one aromatic ring
    aromatic_ring_pattern = Chem.MolFromSmarts("c1ccccc1")
    if not mol.HasSubstructMatch(aromatic_ring_pattern):
        return False, "No aromatic ring found"

    # Check for the presence of at least one oxygen atom in the molecule
    oxygen_atom_pattern = Chem.MolFromSmarts("[O]")
    if not mol.HasSubstructMatch(oxygen_atom_pattern):
        return False, "No oxygen atoms found in the molecule"

    # If all checks pass, the molecule is classified as an anthocyanidin cation
    return True, "Contains flavylium cation structure with oxygenated aromatic rings"