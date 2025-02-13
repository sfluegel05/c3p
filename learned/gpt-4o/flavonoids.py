"""
Classifies: CHEBI:72544 flavonoids
"""
from rdkit import Chem

def is_flavonoids(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    Flavonoids include a wide range of phenyl-substituted benzopyran structures and derivatives.
    
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
        # Flavone basic structure
        Chem.MolFromSmarts("c1cc(O)c2c(c1)C=CC(=O)O2"),    # Flavone
        # Isoflavone structure
        Chem.MolFromSmarts("c1cc2c(cc1)C(=O)C=C(O2)"),     # Isoflavone
        # Chalcone structure
        Chem.MolFromSmarts("c1cc(O)cc(c1)C=CC(=O)c2ccccc2"), # Chalcone
        # Neoflavonoid structure
        Chem.MolFromSmarts("c1ccc2c(c1)C3=CC(=O)C=C(O3)C2"), # Neoflavonoid core
        # Flavanone structure
        Chem.MolFromSmarts("c1ccc2c(c1)C(C2)OC(=O)"),      # Flavanone
        # Flavan-3-ol structure
        Chem.MolFromSmarts("c1cc2c(c(c1)O)C(O)C(C2)O"),    # Flavan-3-ol
        # Anthocyanin structure (flavylium ion)
        Chem.MolFromSmarts("c1ccc2c3c(c(O)c1)C=CC(=O)c3cc(=O)o2"), # Anthocyanidin core
        # Flavonol structure
        Chem.MolFromSmarts("c1cc2c(c(c1)O)C(=O)C=C(O2)"),  # Flavonol
        # Coumestan structure
        Chem.MolFromSmarts("c1ccc2c(c1)C3=CC(=O)OC32"),    # Coumestan
        # Rotenoid structure
        Chem.MolFromSmarts("c1ccc2c3c(c(ccc3O)O2)c1"),     # Rotenoid
    ]
    
    # Check for the presence of any flavonoid core structure
    for pattern in flavonoid_core_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains flavonoid substructure (phenyl-substituted chromanone or related pattern)"
    
    return False, "No core flavonoid substructure found"