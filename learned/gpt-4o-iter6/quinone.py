"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Quinones have a fully conjugated cyclic dione structure, including simple benzoquinones and complex polycyclic forms, derived from aromatic compounds by converting specific groups into carbonyls, covering aromatic and pseudo-aromatic systems.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Expanded SMARTS patterns for quinones, accounting for broader structural motifs
    quinone_patterns = [
        Chem.MolFromSmarts("C1=CC(=O)C=CC(=O)C1"),             # Simple 1,4-benzoquinone
        Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2C1=O"),           # Polycyclic quinone
        Chem.MolFromSmarts("O=C1C(C=CC2=CC=CC=C12)=O"),         # Anthraquinone-like
        Chem.MolFromSmarts("C1=CC2=C(C=C1)C(=O)C=CC2=O"),       # Naphthoquinone-like structures
        Chem.MolFromSmarts("O=C1C=C2C3=CC=CC=C3C2=CC1=O"),      # Tetracyclic-like
        Chem.MolFromSmarts("C1=CC(=O)C=C2C=CC=C3C(=O)OC=C2C3"), # Extended aromatic conjugation
        Chem.MolFromSmarts("O=C1C(=O)C=CC2=C1C=CC(=C2)O"),      # Cross-conjugated
        Chem.MolFromSmarts("O=C1OC=C2C=CC3=CC=COC3C2=C1"),      # Heterocyclic derivatives with oxygen
        Chem.MolFromSmarts("O=C1C(N)=CC=C(O)C1=O"),             # Quinoid nitrogen heterocycle
    ]
    
    # Check for any pattern match
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule contains a fully conjugated cyclic dione structure indicative of a quinone"
    
    return False, "Molecule does not contain a quinone structural feature"