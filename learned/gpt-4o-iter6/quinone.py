"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Quinones have a fully conjugated cyclic dione structure, often derived 
    from aromatic compounds, where an even number of -CH= groups are 
    replaced by -C(=O)- groups within an aromatic ring. Polycyclic and 
    heterocyclic analogues are included.
    
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
    
    # Extended SMARTS patterns for quinones, including polycyclic and heterocyclic variants
    quinone_patterns = [
        Chem.MolFromSmarts("C1=CC(=O)C=CC1=O"),                # Simple 1,4-benzoquinone
        Chem.MolFromSmarts("C1=CC(=O)C=COC1=O"),              # Heterocyclic analogues
        Chem.MolFromSmarts("O=C1C=CC(=O)C=C1"),                # Naphthoquinone-like
        Chem.MolFromSmarts("O=C1C=CC2=C1C=CC(=O)C=C2"),        # Anthraquinone-like
        Chem.MolFromSmarts("O=C1C=CC2=CC=CC=C2C1=O"),          # Extended polycyclic quinone
        Chem.MolFromSmarts("C1(=O)C=CC2=CC(=O)C=CC=C12")       # Another polycyclic variation
    ]
    
    # Check for any pattern match
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule contains a fully conjugated cyclic dione structure indicative of a quinone"
    
    return False, "Molecule does not contain a quinone structural feature"