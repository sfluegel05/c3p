"""
Classifies: CHEBI:35757 monocarboxylic acid anion
"""
"""
Classifies: monocarboxylic acid anion
"""
from rdkit import Chem

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
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the net formal charge of the molecule
    net_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if net_charge != -1:
        return False, f"Molecule has net charge {net_charge}, expected -1"

    # Find all carboxyl groups (protonated or deprotonated)
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=O)[O;H1,-1]")
    carboxyl_groups = mol.GetSubstructMatches(carboxyl_pattern)
    num_carboxyl_groups = len(carboxyl_groups)
    
    if num_carboxyl_groups != 1:
        return False, f"Found {num_carboxyl_groups} carboxyl groups, expected exactly 1"
    
    # Find deprotonated carboxylate groups
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    carboxylate_groups = mol.GetSubstructMatches(carboxylate_pattern)
    num_carboxylate_groups = len(carboxylate_groups)
    
    if num_carboxylate_groups != 1:
        return False, f"Found {num_carboxylate_groups} carboxylate groups, expected exactly 1"

    # Collect indices of carboxylate oxygen atoms
    carboxylate_oxygen_indices = [match[1] for match in carboxylate_groups]

    # Check for other atoms with negative formal charge (excluding carboxylate oxygen)
    negative_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() < 0 and atom.GetIdx() not in carboxylate_oxygen_indices]
    if len(negative_atoms) > 0:
        return False, "Found other negatively charged atoms besides carboxylate oxygen"

    # Check for atoms with positive formal charge
    positive_atoms = [atom for atom in mol.GetAtoms() if atom.GetFormalCharge() > 0]
    if len(positive_atoms) > 0:
        return False, "Molecule contains positively charged atoms (zwitterion)"

    # If all conditions are met, classify as monocarboxylic acid anion
    return True, "Molecule is a monocarboxylic acid anion with exactly one deprotonated carboxyl group and no other ionizable groups"