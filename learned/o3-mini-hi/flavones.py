"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones
A flavone is defined as a flavonoid having a 2-aryl-1-benzopyran-4-one (2-phenylchromen-4-one) skeleton.
This version uses an updated SMARTS pattern to detect the flavone core.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone derivative based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES string into an RDKit molecule.
      2. Generate the Murcko scaffold to remove many peripheral substituents.
      3. Define an updated SMARTS pattern for the 2-phenylchromen-4-one core.
         This pattern represents a phenyl ring attached at the 2-position of a chromen-4-one.
      4. Check for the pattern in both the full molecule and its scaffold.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a flavone derivative, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute the Murcko scaffold (removes peripheral substituents).
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        scaffold = mol  # Fallback on full molecule if scaffold extraction fails.
    
    # Define an updated SMARTS pattern for the 2-phenylchromen-4-one (flavone) core.
    # The pattern "c1ccc(cc1)-c2oc(=O)c3ccccc23" represents:
    #   A phenyl ring (c1ccc(cc1)) attached to a chromen-4-one unit (c2oc(=O)c3ccccc23).
    flavone_pattern = Chem.MolFromSmarts("c1ccc(cc1)-c2oc(=O)c3ccccc23")
    if flavone_pattern is None:
        return False, "Error defining the flavone core SMARTS pattern"
    
    # Check if the full molecule contains the flavone core.
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Molecule contains a 2-phenylchromen-4-one core in its full structure"
    
    # As a fallback, check if the Murcko scaffold contains the flavone core.
    if scaffold.HasSubstructMatch(flavone_pattern):
        return True, "Molecule contains a 2-phenylchromen-4-one core in its Murcko scaffold"
    
    # If neither check finds the expected core, classification fails.
    return False, "Does not contain a 2-phenylchromen-4-one (flavone) core in its structure or scaffold"