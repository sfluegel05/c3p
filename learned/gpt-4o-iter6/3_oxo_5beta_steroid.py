"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is any 3-oxo steroid that has a beta-configuration at position 5.

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

    # Define a SMARTS pattern for the 3-oxo group as part of the steroid structure
    # The pattern considers a carbonyl group (C=O) connected within a cyclopentanoperhydrophenanthrene structure (common core for steroids)
    oxo_steroid_pattern = Chem.MolFromSmarts("[C@@]12([C@]3(CC[C@@H]4[C@]([C@@]3(C1)C)(CCC5[C@@]4(CCC5C(C)C)C)C)C)CCC(=O)[C@]2(C)CCC")
    if not mol.HasSubstructMatch(oxo_steroid_pattern):
        return False, "3-oxo steroid core pattern not found"
    
    # Check for the presence of 5beta stereochemistry
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    # This is a simplification to match chiral centers that relate to position 5 in a relevant steroid context
    # Function may need refining based on specific chirality representation in RDKit for these molecules
    five_beta_stereo = any(idx == 9 and code == 'S' for idx, code in chiral_centers)  # Adjust the idx according to specific position in structure
    
    if not five_beta_stereo:
        return False, "The 5beta stereochemistry was not found or does not match"

    return True, "Molecule is identified as a 3-oxo-5beta-steroid with appropriate stereochemistry"