"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
"""
Classifies: CHEBI:35341 3-oxo-5beta-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid has a ketone at position 3 and beta configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5beta-steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic steroid core (4 fused rings)
    # More flexible pattern that matches the cyclopentanoperhydrophenanthrene core
    steroid_core = Chem.MolFromSmarts("C1C[C@H]2[C@H]3CC[C@H]4CCCC4[C@H]3CC2C1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # Look for C(=O) in first ring
    ketone_pattern = Chem.MolFromSmarts("C1CC(=O)CC2")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group at position 3"

    # Check for 5-beta configuration
    # In 5-beta steroids, the A/B ring junction has the hydrogen in alpha position (down)
    # This creates the characteristic bent shape of 5-beta steroids
    beta_config_pattern = Chem.MolFromSmarts("[C]1[CH2][C](=O)[CH2][C@@H]2")
    if not mol.HasSubstructMatch(beta_config_pattern):
        return False, "Not a 5-beta configuration"

    # Additional validation: Check carbon count (steroids typically have 19-27 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 19 or carbon_count > 30:
        return False, f"Carbon count ({carbon_count}) outside typical steroid range (19-30)"

    # Check for reasonable molecular weight range for steroids
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 250 or mol_weight > 500:
        return False, f"Molecular weight {mol_weight:.1f} outside typical steroid range (250-500)"

    # If we get here, it's a 3-oxo-5beta-steroid
    return True, "Molecule contains steroid core with 3-oxo group and 5-beta configuration"