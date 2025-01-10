"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    A 3-oxo-Delta(4) steroid contains a steroid backbone with a 3-oxo group and a Delta(4) double bond.

    Args:
        smiles (str): SMILES string of the chemical entity.

    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define generalized steroid backbone pattern
    steroid_pattern = Chem.MolFromSmarts('[#6]1[#6][#6][#6]2[#6]([#6]1)[#6][#6][#6]3[#6]2[#6][#6][#6]4[#6]3[#6][#6][#6]4')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Verify the presence of the 3-oxo (C=O) group
    oxo_match = mol.GetSubstructMatches(Chem.MolFromSmarts('C=O'))
    if not oxo_match:
        return False, "No 3-oxo group found"
    
    # Verify the α,β unsaturated bond (Delta(4) bond)
    double_bond_pattern = Chem.MolFromSmarts('[#6]=[#6]')
    if not mol.HasSubstructMatch(double_bond_pattern):
        return False, "No Delta(4) double bond found"
    
    return True, "Molecule classified as 3-oxo-Delta(4) steroid with appropriate moieties"