"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is characterized by an amide group attached to a long carbon chain.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define pattern for amide group
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide linkage found"
    
    # Check for sufficient total number of carbons
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # We now classify a "long" chain as having at least 12 carbon atoms overall
    if total_carbons < 12:
        return False, f"Total carbon count too low for fatty amide, found {total_carbons}"
    
    # Include branched or modified chains by checking for at least some consecutive carbons
    carbon_chain_pattern = Chem.MolFromSmarts("[*]~[*]~[*]~[*]~[*]~[*]~[*]~[*]")  # A more lenient pattern
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "No sufficient consecutive carbon chain detected"

    return True, "Contains an amide linkage and carbon structure indicative of a fatty amide"