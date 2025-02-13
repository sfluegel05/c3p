"""
Classifies: CHEBI:35359 carboxamidine
"""
"""
Classifies: CHEBI:35813 carboxamidine

Carboxamidines are compounds having the structure RC(=NR)NR2. The term is used 
as a suffix in systematic nomenclature to denote the -C(=NH)NH2 group including 
its carbon atom.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxamidine patterns
    carboxamidine_patterns = [
        Chem.MolFromSmarts("[NX3][CX3]([NX3])=[NX2]"),  # Original pattern
        Chem.MolFromSmarts("[NX3][CX3]([NX3])=[NX2]C"),  # Allow additional substitutions
        Chem.MolFromSmarts("[NX3][CX3]([NX3])=[NX2]N"),  # Allow additional substitutions
        # Add more patterns as needed
    ]
    
    for pattern in carboxamidine_patterns:
        if mol.HasSubstructMatch(pattern):
            break
    else:
        return False, "No carboxamidine substructure found"
    
    # Additional checks
    atom_counts = mol.GetAtomWithImplicitHydrogens().GetNumAtomicNumSetOnlyAtoms()
    if atom_counts.get(7, 0) < 3:  # At least 3 nitrogens
        return False, "Insufficient nitrogen atoms for carboxamidine"
    
    if atom_counts.get(6, 0) < 1:  # At least 1 carbon
        return False, "No carbon atom found"
    
    # Check for specific exceptions or edge cases
    # ...
    
    return True, "Molecule contains the carboxamidine substructure R-C(=NR)-NR2"