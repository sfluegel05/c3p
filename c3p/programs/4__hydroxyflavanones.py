"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone has a hydroxy group at the 4' position of the phenyl group B in a flavanone structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a mol object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Define the 4'-hydroxyflavanone pattern generically
    # Allow any B-ring with a 4'-OH group and flexible flavanone structure
    # The pattern represents an aromatic carbon with an OH group in 4' position linked to a flavanone
    flavanone_pattern = Chem.MolFromSmarts("c1ccc(O)cc1[C@H]2CC(=O)c3c(O)cccc3O2")
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No 4'-hydroxyflavanone pattern found"

    return True, "Matches the 4'-hydroxyflavanone structure"

# Example usage
smiles_example = "Oc1ccc(cc1)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1"
is_4__hydroxyflavanones(smiles_example)