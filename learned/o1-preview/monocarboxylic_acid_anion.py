"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion is formed when the carboxy group of a monocarboxylic acid is deprotonated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylate anion group [CX3](=O)[O-]
    carboxylate_anion_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_anion_pattern)
    num_carboxylate = len(carboxylate_matches)

    if num_carboxylate != 1:
        return False, f"Found {num_carboxylate} carboxylate anion groups, need exactly 1"

    # Look for carboxylic acid group [CX3](=O)[OX1H1]
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxylic_acid = len(carboxylic_acid_matches)

    if num_carboxylic_acid > 0:
        return False, f"Found {num_carboxylic_acid} carboxylic acid groups, should have none"

    # Ensure net charge is -1
    charge = Chem.GetFormalCharge(mol)
    if charge != -1:
        return False, f"Molecule has net charge {charge}, should be -1"

    # Check that there are no other acidic functional groups (e.g., additional deprotonated acids)
    # For simplicity, we'll assume monocarboxylic acids only have one anionic group
    # Molecule passes the checks
    return True, "Contains exactly one carboxylate anion group derived from a monocarboxylic acid"