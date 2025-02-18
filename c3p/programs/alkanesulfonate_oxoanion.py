"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
"""
Classifies: CHEBI: alkanesulfonate oxoanion
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is characterized by the presence of at least one sulfonate group (-SO3-)
    attached to a carbon atom, with the carbon potentially being part of a chain or other groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the sulfonate group pattern (-SO3-)
    sulfonate_pattern = Chem.MolFromSmarts("[SX4](=[OX1])(=[OX1])([OX1-])")
    if not mol.HasSubstructMatch(sulfonate_pattern):
        return False, "No sulfonate group (-SO3-) found"

    # Check if the sulfonate group is attached to a carbon atom
    sulfonate_carbon_pattern = Chem.MolFromSmarts("[CX4][SX4](=[OX1])(=[OX1])([OX1-])")
    if not mol.HasSubstructMatch(sulfonate_carbon_pattern):
        return False, "Sulfonate group not attached to a carbon atom"

    # Check for the presence of at least one carbon atom in the molecule
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 1:
        return False, "No carbon atoms found in the molecule"

    # Check for the presence of at least one oxygen atom in the molecule
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen atoms for a sulfonate group"

    # Check for the presence of at least one sulfur atom in the molecule
    s_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
    if s_count < 1:
        return False, "No sulfur atom found in the molecule"

    # Check for the presence of a negative charge on the sulfonate group
    negative_charge_count = sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge() == -1)
    if negative_charge_count < 1:
        return False, "No negative charge found on the sulfonate group"

    # Check for invalid functional groups (e.g., thiols, sulfides)
    invalid_functional_groups = ["[SX2]", "[SX1]"]
    for group in invalid_functional_groups:
        pattern = Chem.MolFromSmarts(group)
        if mol.HasSubstructMatch(pattern):
            return False, f"Molecule contains invalid functional group: {group}"

    return True, "Contains a sulfonate group (-SO3-) attached to a carbon chain or other groups"