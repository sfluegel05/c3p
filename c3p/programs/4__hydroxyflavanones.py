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
    # This pattern represents a flavanone core [C1(COC2=O)] with a 4'-hydroxy group on an aromatic ring linked at C1
    flavanone_pattern = Chem.MolFromSmarts("[OH]c1ccc([C@H]2CC(=O)c3c(cccc3O)O2)cc1")
    
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "No 4'-hydroxyflavanone pattern found"
    
    return True, "Matches the 4'-hydroxyflavanone structure"

# Example usage
smiles_example = "Oc1ccc(cc1)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1"
print(is_4__hydroxyflavanones(smiles_example))