"""
Classifies: CHEBI:83970 cardiac glycoside
"""
from rdkit import Chem

def is_cardiac_glycoside(smiles: str):
    """
    Determines if a molecule is a cardiac glycoside based on its SMILES string.
    Cardiac glycosides typically have a steroid skeleton with a lactone ring
    plus sugar residues.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiac glycoside, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More general steroid backbone definition, accommodating cardiac glycosides
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C3C=CCC4C=CCCC4C3CC2C1')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Detect a range of lactone ring formations
    lactone_patterns = [
        Chem.MolFromSmarts('OC1=CCCCC1'),  # Adaptable lactone structure
        Chem.MolFromSmarts('C1=COC(=O)C1') # Common lactone forms
    ]
    lactone_found = any(mol.HasSubstructMatch(pattern) for pattern in lactone_patterns)
    if not lactone_found:
        return False, "No lactone ring found"
    
    # Extended pattern to identify sugar moieties
    sugar_pattern = Chem.MolFromSmarts('[CX4H1,CX4H2]O[CX4H1,CX4H2]') 
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moieties found"
    
    # Check for more than one sugar-related attachment to ensure complexity
    sugar_attachments = len(mol.GetSubstructMatches(sugar_pattern))
    if sugar_attachments < 2:
        return False, "Insufficient sugar attachments"
    
    return True, "Contains steroid backbone, lactone ring, and glycosidic bonds characteristic of cardiac glycosides."