"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone has a hydroxy group at the 4' position of the phenyl group in a flavanone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the flavanone structure with a 4'-hydroxy group using SMARTS
    # [cH] represents an aromatic carbon with implicit hydrogen (for position flexibility)
    flavanone_pattern = Chem.MolFromSmarts("Oc1ccc(cc1)C2=CC(=O)C3=CC=CC=C3O2")
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No 4'-hydroxyflavanone pattern found"

    return True, "Contains the 4'-hydroxyflavanone core structure"

# Example usage
smiles_example = "Oc1ccc(cc1)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1"
is_4__hydroxyflavanones(smiles_example)