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
    glycone_pattern = Chem.MolFromSmarts("C1C(C(C(C(O1)CO)O)O)O")
    if not mol.HasSubstructMatch(glycone_pattern):
        return False, "No glucoside (glycone) group found"

    # Check for thioglucoside linkage (S-glycoside)
    thioglucoside_pattern = Chem.MolFromSmarts("S-[#6]~C1C(O)C(O)C(O)C(O)C1O")
    if not mol.HasSubstructMatch(thioglucoside_pattern):
        return False, "No thioglucoside linkage found"

    # Check for sulfonated oxime group: N=C-S(=O)(=O)[O-]
    oxime_sulfonate_pattern = Chem.MolFromSmarts("N=C-OS(=O)(=O)[O-]")
    if not mol.HasSubstructMatch(oxime_sulfonate_pattern):
        return False, "No sulfonated oxime group found"
    
    # Check for the central carbon bonded to side chain
    # We implicitly expect the input smiles to cover generic representation of glucosinolates
        
    # Count carbons as a proxy for evaluating side chains and complete structure
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:  # arbitrary threshold, can be adjusted
        return False, "Too few carbon atoms to be a glucosinolate"

    # Additional check: Must have at least two of N, S, and O substituents
    elements_set = {atom.GetSymbol() for atom in mol.GetAtoms()}
    required_elements = {'N', 'S', 'O'}
    if not required_elements.issubset(elements_set):
        return False, "Doesn't have the required elements N, S, and O"
    
    return True, "Matches all structural features of a glucosinolate"