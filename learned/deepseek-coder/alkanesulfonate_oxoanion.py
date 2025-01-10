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
    An alkanesulfonate oxoanion is characterized by the presence of a sulfonate group (-SO3-) 
    attached to a simple carbon chain or other groups.

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

    # Check for the presence of additional functional groups that would disqualify the molecule
    # Exclude molecules with rings, multiple sulfonate groups, or other complex structures
    ring_pattern = Chem.MolFromSmarts("[R]")
    if mol.HasSubstructMatch(ring_pattern):
        return False, "Molecule contains rings, which are not allowed in simple alkanesulfonates"

    # Check for multiple sulfonate groups
    sulfonate_matches = mol.GetSubstructMatches(sulfonate_pattern)
    if len(sulfonate_matches) > 1:
        return False, "Molecule contains multiple sulfonate groups, which is not allowed in simple alkanesulfonates"

    # Check for other functional groups (e.g., amines, alcohols, etc.)
    other_functional_groups = ["[NX3]", "[OX2H]", "[CX3](=O)", "[NX4+]"]
    for group in other_functional_groups:
        pattern = Chem.MolFromSmarts(group)
        if mol.HasSubstructMatch(pattern):
            return False, f"Molecule contains additional functional group: {group}"

    # Check molecular weight to ensure it's not too complex
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 500:
        return False, "Molecular weight too high for a simple alkanesulfonate"

    # Check number of rotatable bonds to ensure it's not too complex
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 10:
        return False, "Too many rotatable bonds for a simple alkanesulfonate"

    return True, "Contains a sulfonate group (-SO3-) attached to a simple carbon chain or other groups"