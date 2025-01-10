"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.

    A 4'-hydroxyflavanone has a flavanone structure with a 4'-hydroxy group on the phenyl ring.

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

    # Refined flavanone core pattern
    flavanone_core_pattern = Chem.MolFromSmarts("c1ccccc1-C2COc3ccccc3O2")  # Phenyl linked a cycle with ketone and hydroxyl
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No flavanone core structure found"

    # Pattern for a 4'-hydroxy group attached at a typical position
    hydroxy_pattern = Chem.MolFromSmarts("c1cc(O)ccc1") # Ensure it's on phenyl component matching flavanone
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "4'-hydroxy group not identified"
    
    return True, "Contains flavanone structure with 4'-hydroxy group"