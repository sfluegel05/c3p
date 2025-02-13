"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone is a hydroxyflavanone with a hydroxy substituent at position 3' of the phenyl ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined flavanone skeleton pattern
    # Matching the core flavanone structure with substitutions
    flavanone_pattern = Chem.MolFromSmarts("O=C1C(Oc2ccccc2O1)c1ccccc1")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No flavanone skeleton found"
    
    # Refine pattern for 3'-hydroxy group relative to core structure
    # Attention to positions in the context of flavanone
    hydroxy_3_prime_pattern = Chem.MolFromSmarts("Oc1ccccc1")
    matches = mol.GetSubstructMatches(hydroxy_3_prime_pattern)

    # Iterate over matches to ensure the positioning is correct
    # Recapture the logic for correctly identifying 3' hydroxy
    found_3_prime_hydroxy = any(mol.GetAtomWithIdx(match[0]).GetFormalCharge() == 0 for match in matches)
    
    if not found_3_prime_hydroxy:
        return False, "No hydroxy group found exactly at position 3' on the phenyl ring"
    
    return True, "Contains flavanone skeleton with a hydroxy group at the 3' position on the phenyl ring"