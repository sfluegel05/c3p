"""
Classifies: CHEBI:24043 flavones
"""
"""
Classifies: Flavones
A flavone is defined as a flavonoid having a 2-aryl-1-benzopyran-4-one (2-phenylchromen-4-one) skeleton.
This version attempts to detect the flavone core by looking both at the Murcko scaffold (to remove 
peripheral substituents) and at the original structure.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_flavones(smiles: str):
    """
    Determines if a molecule is a flavone derivative based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES string into an RDKit molecule.
      2. Compute the Murcko scaffold to (potentially) remove distracting substituents.
      3. Define a SMARTS pattern for a 2-phenylchromen-4-one skeleton.
      4. Check for the pattern in both the full molecule and its scaffold.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a flavone derivative, False otherwise.
        str: Reason for the classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compute the Murcko scaffold for the molecule.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    # If scaffold extraction fails, warn and use the original molecule.
    if scaffold is None:
        scaffold = mol
    
    # Define a SMARTS pattern for the 2-phenylchromen-4-one (flavone) core.
    # It represents a benzopyran-4-one system with an attached phenyl group.
    flavone_pattern = Chem.MolFromSmarts("c1ccccc1-C2=CC(=O)Oc3ccccc23")
    if flavone_pattern is None:
        return False, "Error defining the flavone core SMARTS pattern"
    
    # First, try matching the pattern on the full molecule.
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Molecule contains a 2-phenylchromen-4-one core in its full structure"
    
    # As a fallback, try matching on the Murcko scaffold.
    if scaffold.HasSubstructMatch(flavone_pattern):
        return True, "Molecule contains a 2-phenylchromen-4-one core in its Murcko scaffold"
    
    # If neither search finds the pattern, return False.
    return False, "Does not contain a 2-phenylchromen-4-one (flavone) core in its structure or scaffold"