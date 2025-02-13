"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids typically include a phenyl-substituted benzopyran structure with variations like chromone or chalcone derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is likely a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Extended SMARTS patterns for flavonoid substructures
    flavonoid_core_patterns = [
        # Basic flavonoid structure
        Chem.MolFromSmarts("c1cc(O)c2c(c1)C=CC(=O)O2"),    # Flavone
        Chem.MolFromSmarts("c1cc2c(cc1)C(=O)C=C(O2)"),     # Isoflavone
        # Chalcone structure
        Chem.MolFromSmarts("c1cc(O)cc(c1)C=CC(=O)c2ccccc2"), # Chalcone
        # Neoflavonoid structure
        Chem.MolFromSmarts("c1ccccc1C2=CC(=O)C3=CC=CC(O)=C23"), # Neoflavonoid
        # Flavanone structure
        Chem.MolFromSmarts("c1ccc2c(c1)C(C2)OC(=O)"),      # Flavanone
        # Additional pattern to cover more classes
        Chem.MolFromSmarts("c1cc2c(c(c1)O)C(=C(C(=O)c3ccc(O)cc3)O2)"), # Flavan-3-ol
    ]
    
    # Check for the presence of any flavonoid core structure
    for pattern in flavonoid_core_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains flavonoid substructure (phenyl-substituted chromanone or related pattern)"
    
    return False, "No core flavonoid substructure found"