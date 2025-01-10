"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_1_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(1) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for the steroid backbone with a 3-keto group and a Delta(1) (1-position) double bond
    steroid_pattern = Chem.MolFromSmarts("C=1C[C@H]2[C@@H](C(C)=O)CCC3=CCC=C[C@]3(C)[C@@H]2CC1")
    
    # Check for the steroid pattern in the molecule
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "3-oxo-Delta(1) steroid structure not found"

    return True, "Molecule contains a 3-oxo-Delta(1) steroid structure"