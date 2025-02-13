"""
Classifies: CHEBI:38757 isoflavones
"""
"""
Classifies: CHEBI:28845 isoflavone

Any isoflavonoid with a 3-aryl-1-benzopyran-4-one (3-aryl-4H-chromen-4-one) skeleton 
and its substituted derivatives.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_isoflavone(smiles: str):
    """
    Determines if a molecule is an isoflavone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the 3-aryl-4H-chromen-4-one skeleton
    isoflavone_pattern = Chem.MolFromSmarts("[o,2]1[#6]:2:[#6]:[#6]:[#6]:[#6]:2:[#6]:1-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1")
    if not mol.HasSubstructMatch(isoflavone_pattern):
        return False, "Does not contain the 3-aryl-4H-chromen-4-one skeleton"
    
    # Check for aromaticity
    if not Chem.AromaticityCalculator.AromaticAtomsCount(mol, 6):
        return False, "Molecule is not aromatic"
    
    # Check for presence of carbonyl group
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "Does not contain a carbonyl group"
    
    # Check for presence of heterocyclic oxygen
    hetero_oxygen_pattern = Chem.MolFromSmarts("[O;r]")
    if not mol.HasSubstructMatch(hetero_oxygen_pattern):
        return False, "Does not contain a heterocyclic oxygen"
    
    return True, "Molecule contains the 3-aryl-4H-chromen-4-one skeleton and exhibits features of an isoflavone"