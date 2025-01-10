"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone has a flavanone core structure with a hydroxy group
    on the phenyl ring specifically at the 3' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the flavanone core structure involving chroman-4-one moiety
    flavanone_core_pattern = Chem.MolFromSmarts("O=C1CC(c2ccccc2)Oc2ccccc12")
    if not mol.HasSubstructMatch(flavanone_core_pattern):
        return False, "Flavanone core structure not found"

    # Define SMARTS for 3'-position hydroxy specifically on the A ring, 
    # looking for hydroxy on the meta position of the phenyl group
    hydroxy_3_prime_pattern = Chem.MolFromSmarts("c1cc(O)ccc1")
    
    # Find all matching phenyl groups
    phenyl_ring = Chem.MolFromSmarts("c1ccccc1")
    phenyl_matches = mol.GetSubstructMatches(phenyl_ring)
    found_3_prime_hydroxy = False

    for match in phenyl_matches:
        submol = Chem.PathToSubmol(mol, match)
        if submol.HasSubstructMatch(hydroxy_3_prime_pattern):
            found_3_prime_hydroxy = True
            break

    if not found_3_prime_hydroxy:
        return False, "No hydroxy group at 3' position on the phenyl ring"

    return True, "Contains flavanone core with a hydroxy group at 3' position"