"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:25438 nucleobase analogue
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    A nucleobase analogue is a molecule that can substitute for a normal nucleobase in nucleic acids.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a nucleobase analogue, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for purine and pyrimidine rings
    purine_pattern = Chem.MolFromSmarts("n1cnc2c(ncnc12)")  # Basic purine ring
    pyrimidine_pattern = Chem.MolFromSmarts("c1cncnc1")      # Basic pyrimidine ring

    # Check for purine ring system
    purine_match = mol.HasSubstructMatch(purine_pattern)
    
    # Check for pyrimidine ring system
    pyrimidine_match = mol.HasSubstructMatch(pyrimidine_pattern)
    
    # Additional check for modified nucleobases (e.g., substitutions or analogues)
    # We can look for heterocyclic rings with nitrogen atoms
    heterocycle_with_N = Chem.MolFromSmarts("[n;R]")  # Any nitrogen in a ring

    heterocycle_N_matches = mol.GetSubstructMatches(heterocycle_with_N)
    has_heterocycle_N = len(heterocycle_N_matches) > 0

    # Check if molecule is aromatic
    is_aromatic = any(atom.GetIsAromatic() for atom in mol.GetAtoms())

    # Decide classification based on matches
    if purine_match or pyrimidine_match:
        return True, "Contains purine or pyrimidine ring system"
    elif has_heterocycle_N and is_aromatic:
        return True, "Contains aromatic heterocycle(s) with nitrogen, characteristic of nucleobase analogues"
    else:
        return False, "Does not contain structures characteristic of nucleobase analogues"