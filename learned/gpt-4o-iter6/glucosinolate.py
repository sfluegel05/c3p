"""
Classifies: CHEBI:24279 glucosinolate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glucosinolate(smiles: str):
    """
    Determines if a molecule is a glucosinolate based on its SMILES string.
    Discriminative features:
    - Contains a glycone group, likely glucose.
    - A thioglucoside linkage (S-glycoside).
    - A central carbon linked via S and N, with an oxime and a sulfonate group.
    - A side chain on the central carbon.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a glucosinolate, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the glycone group using a glucose-like pattern
    # The pattern is now more relaxed to include thioglycosides.
    glycone_pattern = Chem.MolFromSmarts("C1C(C(C(O1)CO)O)O")
    if not mol.HasSubstructMatch(glycone_pattern):
        return False, "No glucoside (glycone) group found"

    # Check for thioglucoside linkage (S-glycoside)
    # This pattern is simplified to look for any sugar linked with a sulfur.
    thioglucoside_pattern = Chem.MolFromSmarts("S-[#6]")
    if not mol.HasSubstructMatch(thioglucoside_pattern):
        return False, "No thioglucoside linkage found"

    # Check for sulfonated oxime group: N=C-OS(=O)(=O)[O-]
    oxime_sulfonate_pattern = Chem.MolFromSmarts("N=C-OS(=O)(=O)[O-]")
    if not mol.HasSubstructMatch(oxime_sulfonate_pattern):
        return False, "No sulfonated oxime group found"
    
    # Check for the central carbon bonded to side chain
    central_c_pattern = Chem.MolFromSmarts("[#6]-[S]-[C]")
    if not mol.HasSubstructMatch(central_c_pattern):
        return False, "Central C structure not well-configured for glucosinolate"

    # Count the elements, ensuring appropriate diversity
    elements_set = {atom.GetSymbol() for atom in mol.GetAtoms()}
    required_elements = {'N', 'S', 'O'}
    if not required_elements.issubset(elements_set):
        return False, "Doesn't have the required elements N, S, and O"
    
    return True, "Matches all structural features of a glucosinolate"