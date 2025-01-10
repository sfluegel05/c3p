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

    # More general steroid core pattern - four fused rings
    # This pattern represents the basic connectivity without being too strict about stereochemistry
    steroid_core = Chem.MolFromSmarts("C1CC2CCC3C(CCC4CCCC43C)C2C1")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # More specific pattern that ensures the ketone is at position 3
    # The pattern includes the A ring and part of the B ring to ensure correct position
    ketone_pattern = Chem.MolFromSmarts("C1CC(=O)CC2CCC3")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group at position 3"

    # Check for 5-beta configuration
    # In 5-beta steroids, the A/B ring junction has the hydrogen in alpha position
    # This creates a cis fusion between rings A and B
    # The pattern looks for the specific stereochemistry at position 5
    beta_config_pattern = Chem.MolFromSmarts("[C]1[CH2][C](=O)[CH2][C@@H]2[CH2][CH2]")
    
    # Alternative 5-beta pattern that might catch some other representations
    beta_config_pattern2 = Chem.MolFromSmarts("[C]1[CH2][C](=O)[CH2][C@]2([H])")
    
    if not (mol.HasSubstructMatch(beta_config_pattern) or mol.HasSubstructMatch(beta_config_pattern2)):
        return False, "Not a 5-beta configuration"

    # Additional validation: Check carbon count (steroids typically have 19-27 carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 19 or carbon_count > 30:
        return False, f"Carbon count ({carbon_count}) outside typical steroid range (19-30)"

    # Check for reasonable molecular weight range for steroids
    mol_weight = Chem.Descriptors.ExactMolWt(mol)
    if mol_weight < 250 or mol_weight > 500:
        return False, f"Molecular weight {mol_weight:.1f} outside typical steroid range (250-500)"

    # If all checks pass, we have a 3-oxo-5beta-steroid
    return True, "Molecule contains steroid core with 3-oxo group and 5-beta configuration"