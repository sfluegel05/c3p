"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    Catechins have a flavan-3-ol skeleton and may have various hydroxyl and other substitutions.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more comprehensive SMARTS pattern for the flavan-3-ol skeleton
    # Flavan-3-ol basic structure: two benzene rings connected by a pyran ring with hydroxyl groups
    flavan3ol_pattern = Chem.MolFromSmarts("C1C[C@H](Oc2c(cccc2O)C1O)Oc3cc(O)cc(O)c3") 
    if not mol.HasSubstructMatch(flavan3ol_pattern):
        return False, "No flavan-3-ol skeleton found"
    
    # Check for presence of hydroxyl groups
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if hydroxyl_count < 2:
        return False, f"Too few hydroxyl groups for a catechin, found {hydroxyl_count}"

    return True, "Molecule matches flavan-3-ol skeleton and has sufficient hydroxyl groups"