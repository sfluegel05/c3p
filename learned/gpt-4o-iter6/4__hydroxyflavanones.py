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
    
    # Define a generic flavanone core pattern
    flavanone_core_pattern = Chem.MolFromSmarts("O=C1CC(Oc2ccccc2)C(O1)c3ccccc3")
    
    # Check for the presence of the flavanone core structure
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "No flavanone core structure identified"
    
    # Define pattern for 4'-hydroxy group on an aromatic ring adjacent to the flavanone core
    four_prime_hydroxy_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    
    # Check for a 4'-OH group on the phenyl ring
    atom_indices = mol.GetSubstructMatch(four_prime_hydroxy_pattern)
    if atom_indices:
        return True, "Contains 4'-hydroxy group on aromatic ring of flavanone core"
    
    return False, "4'-hydroxy group not found at the correct position"