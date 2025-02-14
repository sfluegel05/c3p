"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monocarboxylic_acid_anion(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a monocarboxylic acid anion based on its SMILES string.
    A monocarboxylic acid anion contains a single deprotonated carboxyl group (C(=O)[O-]).

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

    # Count carboxylate groups
    carboxylate_pattern = Chem.MolFromSmarts("C(=O)[O-]")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    # Count carboxylic acid groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O[H]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    

    if len(carboxylate_matches) != 1:
      return False, f"Molecule has {len(carboxylate_matches)} carboxylate groups, it should have exactly 1"

    if len(carboxylic_acid_matches) > 0:
      return False, "Molecule has protonated carboxylic acid groups"
    
    charge = rdMolDescriptors.CalcFormalCharge(mol)

    if charge != -1:
        return False, f"Molecule has a charge of {charge}, should have a charge of -1"

    return True, "Molecule is a monocarboxylic acid anion"