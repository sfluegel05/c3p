"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone has a hydroxy group at the 4' position on the B-ring of a flavanone core structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Define a refined SMARTS pattern for 4'-hydroxyflavanone
    # The flavanone core pattern with hydroxy on the B-ring is targeted
    flavanone_core = Chem.MolFromSmarts("C1(C(C=O)c2ccccc2O1)c3ccc(O)cc3")
    
    # Check the stereochemistry around the flavanone core
    flavanone_stereo_pattern = Chem.MolFromSmarts("[C@H]1(C(=O)c2c(cccc2O1)c3ccc(O)cc3)C")

    # Ensure it has both flavanone core and the specific stereochemistry
    if not mol.HasSubstructMatch(flavanone_core):
        return False, "No flavanone core structure found"
    
    if not mol.HasSubstructMatch(flavanone_stereo_pattern):
        return False, "Incorrect stereochemistry for 4'-hydroxyflavanone"

    return True, "Matches the 4'-hydroxyflavanone structure"

# Example usage
smiles_example = "Oc1ccc(cc1)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1"
print(is_4__hydroxyflavanones(smiles_example))