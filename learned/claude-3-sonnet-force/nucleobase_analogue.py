"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: CHEBI:37585 nucleobase analogue
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
    
    # Define SMARTS patterns for common nucleobase analogue structural features
    purine_pattern = Chem.MolFromSmarts("c1nc2[nH]cnc2[nH]c1")
    pyrimidine_pattern = Chem.MolFromSmarts("c1c[nH]c(=O)[nH]c1")
    amino_pattern = Chem.MolFromSmarts("[NH2,NH]")
    oxo_pattern = Chem.MolFromSmarts("C=O")
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    
    # Check for purine or pyrimidine ring
    if not mol.HasSubstructMatch(purine_pattern) and not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No purine or pyrimidine ring found"
    
    # Check for amino, oxo, and hydroxyl functional groups
    has_amino = mol.HasSubstructMatch(amino_pattern)
    has_oxo = mol.HasSubstructMatch(oxo_pattern)
    has_hydroxyl = mol.HasSubstructMatch(hydroxyl_pattern)
    if not (has_amino or has_oxo or has_hydroxyl):
        return False, "Missing amino, oxo, or hydroxyl functional groups"
    
    # Check for typical size range of nucleobase analogues
    n_heavy_atoms = mol.GetNumHeavyAtoms()
    if n_heavy_atoms < 5 or n_heavy_atoms > 25:
        return False, "Size deviates too much from typical nucleobases"
    
    # Check for aromaticity and planarity
    if not mol.GetNumBEdges() or not AllChem.PlanarityAtomsSelfIter(mol):
        return False, "Molecule is not sufficiently aromatic or planar"
    
    # If all conditions are met, classify as nucleobase analogue
    return True, "Contains purine or pyrimidine ring and typical functional groups of nucleobase analogues"