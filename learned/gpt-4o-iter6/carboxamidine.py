"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    Carboxamidine has the structure RC(=NR)NR2, where R can be hydrogen or other substituents.

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

    # SMARTS patterns for carboxamidine: RC(=NR)NR2, accommodating variable R groups
    # Pattern for various resonance forms and bonding configurations of carboxamidine
    carboxamidine_pattern1 = Chem.MolFromSmarts("[CX3](=[NX2])[NX3H2,NX3H][#6,#1,#7]")  # Covers more bonds to nitrogen
    carboxamidine_pattern2 = Chem.MolFromSmarts("[N]=[C]-[N]")                          # Covers C=N=N type bonding
    carboxamidine_pattern3 = Chem.MolFromSmarts("[CX3]=[NX2][NX3H1]")                   # Other tautomer forms

    # Check the molecule for matching any of the carboxamidine patterns
    if (mol.HasSubstructMatch(carboxamidine_pattern1) or
        mol.HasSubstructMatch(carboxamidine_pattern2) or
        mol.HasSubstructMatch(carboxamidine_pattern3)):
            return True, "Contains a carboxamidine structure: RC(=NR)NR2"
    
    return False, "Does not contain a carboxamidine structure"