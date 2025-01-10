"""
Classifies: CHEBI:36249 bile acid conjugate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate involves a bile acid structure attached to hydrophilic/charged 
    groups such as glycine, taurine, sulfate, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generalized steroid backbone pattern
    steroid_backbone_pattern = Chem.MolFromSmarts("C1CCC2C3C4CC(CCC4)C(O)C3CCC2C1")  
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No bile acid core structure found"

    # Expanded conjugate patterns
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)[O,*,1]")  # Allow variable attachment
    taurine_pattern = Chem.MolFromSmarts("S(=O)(=O)[A,N][*,1]C(=O)")  # Sulfonamide link
    sulfate_pattern = Chem.MolFromSmarts("OS(=O)(=O)[A,*,1]")  # Ester or charged
    glucuronic_acid_pattern = Chem.MolFromSmarts("OC(CO)C(O)C(=O)[O,A][*,1]")  # Variable attachment

    # Match any of the conjugal patterns
    has_glycine = mol.HasSubstructMatch(glycine_pattern)
    has_taurine = mol.HasSubstructMatch(taurine_pattern)
    has_sulfate = mol.HasSubstructMatch(sulfate_pattern)
    has_glucuronic_acid = mol.HasSubstructMatch(glucuronic_acid_pattern)

    if not any([has_glycine, has_taurine, has_sulfate, has_glucuronic_acid]):
        return False, "No known bile acid conjugate pattern found"

    # If both core and conjugate are present
    return True, "Matches bile acid core with conjugate pattern"

# Example SMILES to test
example_smiles = "C1[C@@]2([C@@H](C[C@H]1O)CC[C@H]1[C@@]2(C)[C@@H]2CC[C@@H]3[C@@H](C(CC[C@]23C)=O)NCCS(=O)(=O)O)[C@H](C)O"
result, reason = is_bile_acid_conjugate(example_smiles)
print(result, reason)