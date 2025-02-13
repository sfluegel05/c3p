"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: 3'-hydroxyflavanones
Compounds with a flavanone core structure and a hydroxy substituent at position 3' of the phenyl ring
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # First check for basic flavanone core (chroman-4-one)
    # Note: The =O must be explicit to avoid matching flavones
    flavanone_core = Chem.MolFromSmarts("[#6]1-[#6]-C(=O)-c2c(O1)cccc2")
    if flavanone_core is None:
        return False, "Error in SMARTS pattern"
    
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "Missing flavanone core structure"

    # Check for phenyl ring at position 2 with hydroxyl at position 3'
    # The phenyl ring must be connected to the carbon between the O and C=O
    # [#6H1,#6H0] allows for substituted carbons
    pattern = Chem.MolFromSmarts("""
        [#6]1-[#6]-C(=O)-c2c(O1)cccc2
        [$([#6H1]),$([#6H0])]1@[#6](-[#6]-C(=O)-[#6]2-[#6]-[#6]-[#6]-[#6]-[#6]-2-O-1)
        =[$([#6H1]),$([#6H0])]-[#6](-[OH1])=[$([#6H1]),$([#6H0])]-[#6]=[$([#6H1]),$([#6H0])]
    """)
    
    if pattern is None:
        return False, "Error in SMARTS pattern"

    if mol.HasSubstructMatch(pattern):
        return True, "Contains flavanone core with 3'-hydroxy substituent"
    
    # If we get here, molecule has flavanone core but no 3'-OH
    return False, "Has flavanone core but missing 3'-hydroxy substituent"