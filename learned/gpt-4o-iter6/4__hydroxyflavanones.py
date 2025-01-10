"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent located at position 4' on the phenyl group.

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
    
    # Define a more general flavanone core pattern
    # Flavanone core consists of a 2-phenylchroman-4-one structure
    flavanone_core_pattern = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc3ccccc13")
    
    # Check for the presence of the flavanone core structure
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No flavanone core structure identified"
    
    # Define pattern for 4'-hydroxy group on an aromatic ring attached to the flavanone core
    # Ensures OH group is para to the attachment point
    four_prime_hydroxy_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    
    # Check for a 4'-OH group on the phenyl ring
    atom_indices = mol.GetSubstructMatch(four_prime_hydroxy_pattern)
    if atom_indices:
        # Verify that the OH group is on the correct phenyl ring
        phenyl_ring = Chem.MolFromSmarts("c1ccccc1")
        phenyl_match = mol.GetSubstructMatch(phenyl_ring)
        if phenyl_match and any(idx in atom_indices for idx in phenyl_match):
            return True, "Contains 4'-hydroxy group on aromatic ring of flavanone core"
    
    return False, "4'-hydroxy group not found at the correct position"