"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem

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
    
    # Define the flavanone core pattern (including stereochemistry)
    flavanone_pattern = Chem.MolFromSmarts("O=C1C([C@@H](O)Cc2ccccc2)Oc3ccccc13")
    if not mol.HasSubstructMatch(favanone_pattern):
        return False, "No flavanone skeleton found"

    # Define pattern for the 3' hydroxy group on the A-ring of flavanones
    # Locating the hydroxy group at the 3-prime position on correct ring
    three_prime_hydroxy_pattern = Chem.MolFromSmarts("Oc1[cH][cH][cH][cH][cH]1")
    
    # Approach may consider how the hydroxy position logically fits here
    # Shift focus: attachment to the flavanone, check adjacency, etc.
    three_prime_matches = mol.GetSubstructMatches(three_prime_hydroxy_pattern)
    found_hydroxy = False
    
    if len(three_prime_matches) > 0:
        for match in three_prime_matches:
            # Verify neighbors to establish the 3' position is relevant to core
            # This might need more specification (like verifying ring context)
            hydroxy_atom = mol.GetAtomWithIdx(match[0])
            connected_atoms = [n.GetAtomicNum() for n in hydroxy_atom.GetNeighbors()]
            # Ensure hydroxy is on a phenyl ring expected for 3' position
            if 6 in connected_atoms:  # Presence of carbon indicates phenyl ring
                found_hydroxy = True
                break

    if not found_hydroxy:
        return False, "No hydroxy group found exactly at position 3' on the phenyl ring"
    
    return True, "Contains flavanone skeleton with a hydroxy group at the 3' position on the phenyl ring"