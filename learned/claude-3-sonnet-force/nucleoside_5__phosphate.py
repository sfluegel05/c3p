"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:37563 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate is a ribosyl or deoxyribosyl derivative of a pyrimidine or purine base
    in which C-5 of the ribose ring is mono-, di-, tri- or tetra-phosphorylated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for pyrimidine or purine base
    purine_pattern = Chem.MolFromSmarts("[*]1[*]c2[nH]c[nH]c2[*]c1[*]")  # Purine
    pyrimidine_pattern = Chem.MolFromSmarts("[*]1[*]c[nH]c([nH]c1[*])[*]")  # Pyrimidine
    base_pattern = purine_pattern | pyrimidine_pattern
    if not mol.HasSubstructMatch(base_pattern):
        return False, "No purine or pyrimidine base found"
    
    # ... (rest of the code remains the same) ...