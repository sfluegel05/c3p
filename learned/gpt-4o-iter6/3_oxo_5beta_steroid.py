"""
Classifies: CHEBI:1624 3-oxo-5beta-steroid
"""
from rdkit import Chem

def is_3_oxo_5beta_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5beta-steroid based on its SMILES string.
    A 3-oxo-5beta-steroid is any 3-oxo steroid that has beta-configuration at position 5.

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
    
    # Define more precise SMARTS for the 3-oxo group within a steroid backbone
    oxo_pattern = Chem.MolFromSmarts("C1(C(=O))[C]([C]2)[C]([C]CCC3)[C](C4)[C]([C]([C](=C4)CC3)C2)=C1")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No 3-oxo group found in the expected steroid structure"

    # General pattern for steroid backbones (adapt for generality)
    steroid_core_pattern = Chem.MolFromSmarts("C1(CCC2C3(C4CCCCC4)CCC3CCC2C1)")
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "Steroid core not found"

    # Check 5beta stereochemistry
    # Get stereocenter indices and evaluate CIP configuration at these indices
    chiral_centers = Chem.FindMolChiralCenters(mol, force=True, includeUnassigned=True)
    stereochemically_correct = any(code == 'S' and mol.GetAtomWithIdx(idx).GetAtomicNum() == 6
                                   for idx, code in chiral_centers if idx == 5)
    
    if not stereochemically_correct:
        return False, "5beta stereochemistry not resolved"

    return True, "Molecule is identified as a 3-oxo-5beta-steroid with appropriate stereochemistry"