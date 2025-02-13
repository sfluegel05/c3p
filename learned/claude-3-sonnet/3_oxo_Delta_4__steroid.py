"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
"""
Classifies: CHEBI:38034 3-oxo-Delta(4) steroid

A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for steroid backbone pattern
    # [C@@]1(CC[C@H]2[C@H]3CCC4=CC(=O)C=C[C@]4(C)[C@H]3C[C@]12C)[H]
    steroid_pattern = Chem.MolFromSmarts("[C@@]1(CC[C@H]2[C@H]3CCC4=CC(=O)C=C[C@]4(C)[C@H]3C[C@]12C)[H]")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Look for 3-oxo group
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C@@]")
    if len(mol.GetSubstructMatches(oxo_pattern)) != 1:
        return False, "No 3-oxo group found or more than one present"
    
    # Look for alpha,beta conjugated C=C bond
    conjugated_pattern = Chem.MolFromSmarts("C=C[C@@]C(=O)")
    if not mol.HasSubstructMatch(conjugated_pattern):
        return False, "No conjugated alpha,beta C=C bond found"

    return True, "Contains 3-oxo group and conjugated alpha,beta C=C bond on steroid backbone"