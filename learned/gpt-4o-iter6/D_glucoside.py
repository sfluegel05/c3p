"""
Classifies: CHEBI:35436 D-glucoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_D_glucoside(smiles: str):
    """
    Determines if a molecule is a D-glucoside based on its SMILES string.
    A D-glucoside has a glucose moiety with the specific stereochemistry
    corresponding to D-glucose.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a D-glucoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glucose moiety pattern with expected stereochemistry for D-glucose
    glucose_pattern = Chem.MolFromSmarts("[C@@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H](CO)O1")
    if not mol.HasSubstructMatch(glucose_pattern):
        return False, "No D-glucopyranose unit found"
    
    # Additional validation can be added to confirm the presence of 
    # a glycosidic linkage, depending on the domain knowledge available
    
    return True, "Contains a D-glucose moiety with proper stereochemistry"

# Example usage
print(is_D_glucoside("CC(C)=CCc1cc(cc(O)c1O)-c1oc2cc(O[C@@H]3O[C@H](CO)[C@@H](O)[C@H](O)[C@H]3O)cc(O)c2c(=O)c1CO")) 
# Output should be: (True, 'Contains a D-glucose moiety with proper stereochemistry')