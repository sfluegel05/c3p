"""
Classifies: CHEBI:20156 3-oxo-Delta(1) steroid
"""
"""
Classifies: CHEBI:37663 3-oxo-Delta(1) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_Delta_1__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(1) steroid based on its SMILES string.
    A 3-oxo-Delta(1) steroid is defined as any 3-oxo steroid that contains a double bond between positions 1 and 2.

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

    # Check for steroid backbone (4 fused rings)
    steroid_pattern = Chem.MolFromSmarts("[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)C=C[C@]4(C)[C@H]3CC[C@]12C")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 3-oxo group (carbonyl at position 3)
    oxo_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not any(match[0] == 2 for match in oxo_matches):  # Assuming position 3 is atom index 2
        return False, "No 3-oxo group found"

    # Check for Delta(1) double bond (double bond between positions 1 and 2)
    delta1_pattern = Chem.MolFromSmarts("[C]=[C]")
    delta1_matches = mol.GetSubstructMatches(delta1_pattern)
    if not any(match[0] == 0 and match[1] == 1 for match in delta1_matches):  # Assuming positions 1 and 2 are atom indices 0 and 1
        return False, "No Delta(1) double bond found"

    return True, "Contains steroid backbone with 3-oxo group and Delta(1) double bond"