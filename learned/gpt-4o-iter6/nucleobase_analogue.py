"""
Classifies: CHEBI:67142 nucleobase analogue
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Purine and pyrimidine SMARTS patterns
    purine_pattern = Chem.MolFromSmarts("n1cnc2[nH]cnc12")
    pyrimidine_pattern = Chem.MolFromSmarts("c1ncnc1")
    
    # Check for core nucleobase structures
    has_purine_structure = mol.HasSubstructMatch(purine_pattern)
    has_pyrimidine_structure = mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not (has_purine_structure or has_pyrimidine_structure):
        return False, "No purine or pyrimidine backbone detected"
    
    # Check for typical functional groups or modifications
    # Example modifications: halogens, ketones, amines, hydroxyl, etc.
    nucleobase_modifications = Chem.MolFromSmarts("[OH0,Cl,Br,I,F,NH1,NH2,=O,S]") # rough example pattern for modifications
    if mol.HasSubstructMatch(nucleobase_modifications):
        return True, "Core nucleobase structure present with known modifications"

    return True, "Core nucleobase structure without major natural deviations"

# Examples provided can be tested using this function
# is_nucleobase_analogue('SMILES_STRING_HERE')