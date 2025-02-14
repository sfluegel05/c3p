"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: CHEBI:26266 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is an organic aromatic compound with a structure based on a phenylpropane skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phenylpropane backbone (benzene ring with alkyl chain of 1-6 carbons)
    backbone_pattern = Chem.MolFromSmarts("[c]1ccccc1C[C@@](CC)(CC)CC")
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No phenylpropane backbone found"
    
    # Check for common substituents (hydroxy, methoxy, ester, prenyl, etc.)
    substituents = ['O', 'OC', 'C(=O)O', 'C(=O)OC', 'c1ccccc1', 'CC=C(C)C']
    substituent_pattern = Chem.MolFromSmarts("|".join(["[$({})]/[#6]".format(sub) for sub in substituents]))
    if not mol.HasSubstructMatch(substituent_pattern):
        return False, "No common phenylpropanoid substituents found"
    
    # Look for common phenylpropanoid scaffolds (flavonoids, coumarins, lignans, etc.)
    scaffold_patterns = [
        Chem.MolFromSmarts("c1cc2c(cc1O)OC(C=3C(=O)OC(c3=O)c4ccccc4)CC2"),  # flavonoid
        Chem.MolFromSmarts("C1=C(C(=O)OC1)c2ccccc2"),  # coumarin
        Chem.MolFromSmarts("c1ccc(C[C@H](C)Cc2ccc(OC)cc2)cc1OC")  # lignan
    ]
    if any(mol.HasSubstructMatch(pattern) for pattern in scaffold_patterns):
        return True, "Contains common phenylpropanoid scaffold"
    
    return True, "Contains phenylpropane backbone with aromatic ring and common substituents"