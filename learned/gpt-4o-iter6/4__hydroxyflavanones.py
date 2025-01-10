"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent located at position 4'.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved SMARTS pattern for flavanone core structure
    # Including the heterocyclic core and connection to a phenyl group
    flavanone_core_pattern = Chem.MolFromSmarts("O=C1CC(=O)c2ccccc2O1")  # Core flavanone structure without stereochemistry
    
    # Check for flavanone core in molecule
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No flavanone core structure identified"

    # Detect 4'-hydroxy group, attached to the phenyl ring bonded to flavanone
    four_prime_hydroxy_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    
    # Check if the phenyl ring containing the 4' hydroxy group is connected to the flavanone core
    for match in mol.GetSubstructMatches(four_prime_hydroxy_pattern):
        atom_indices = set(match)
        for _, item in enumerate(mol.GetSubstructMatches(flavanone_core_pattern)):
            flavanone_atoms = set(item)
            if atom_indices & flavanone_atoms:
                return True, "Contains 4'-hydroxyflavanone core and a 4'-hydroxy group"

    return False, "4'-hydroxy group not found at the correct position"