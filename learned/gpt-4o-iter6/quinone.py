"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Quinones have a fully conjugated cyclic dione structure, typically derived 
    from aromatic compounds. We include polycyclic and heterocyclic variants, 
    marked by the presence of two adjacent carbonyl groups within an extended 
    conjugated system (aromatic or pseudo-aromatic).
    
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
    
    # Expanded SMARTS patterns for quinones, including variations and analogues
    quinone_patterns = [
        Chem.MolFromSmarts("C1=CC(=O)C=CC1=O"),                # Simple 1,4-benzoquinone
        Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2C1=O"),          # Extended polycyclic quinone
        Chem.MolFromSmarts("O=C1C(C=CC2=CC=CC=C12)=O"),        # Anthraquinone-like,
        Chem.MolFromSmarts("C1=CC2=C(C=C1)C(=O)C=CC2=O"),      # Naphthoquinone-like structure
        Chem.MolFromSmarts("C1=CC(=O)C=C2C=CC=CC2=C1"),        # Polycyclic with substitution
        Chem.MolFromSmarts("C1=COC(=O)C=C1"),                  # Heterocyclic derivatives
        Chem.MolFromSmarts("C1=CC(=O)C=CC(=O)C1"),             # Extended aromatic with two carbonyls
        Chem.MolFromSmarts("O=C1C=C2C(=O)c3ccccc3C2=CC1"),     # Similar to tetracenequinone structures
    ]
    
    # Check for any pattern match
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule contains a fully conjugated cyclic dione structure indicative of a quinone"
    
    return False, "Molecule does not contain a quinone structural feature"