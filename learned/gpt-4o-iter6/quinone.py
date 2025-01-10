"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Quinones have a fully conjugated cyclic dione structure, typically derived 
    from aromatic compounds, where an even number of -CH= groups are 
    replaced by -C(=O)- groups within an aromatic ring. This includes polycyclic
    and heterocyclic analogues.
    
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
    
    # Extended SMARTS patterns for quinone structures, including polycyclic variants
    # Basic benzoquinone-like pattern
    quinone_patterns = [
        Chem.MolFromSmarts("C1=CC(=O)C=C(C=1)=O"),  # Simple 1,4-benzoquinone
        Chem.MolFromSmarts("C1=CC(=O)C=CC=C1(=O)"), # Extended to allow extra conjugation paths
        Chem.MolFromSmarts("O=C1C=CC2=CC(O)=CC=C2C1"), # Naphthoquinone-like
        Chem.MolFromSmarts("O=C1C=CC2=C1C(=O)C=CC=C2")  # Anthraquinone-like
    ]
    
    # Check for any pattern match
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule contains a fully conjugated cyclic dione structure indicative of a quinone"
    
    return False, "Molecule does not contain a quinone structural feature"