"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    Quinones have a fully conjugated cyclic dione structure, which can be simple benzoquinones or complex polycyclic forms derived from aromatics by converting specific groups into carbonyls, covering aromatic and pseudo-aromatic systems.
    
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

    # SMARTS patterns for quinones
    quinone_patterns = [
        "[O]=C1C=CC(=O)C=C1",                # Simple 1,4-benzoquinone
        "[O]=C1C2=CC=CC=C2C(=O)C=C1",        # Naphthoquinone-like
        "[O]=C1CC=CC(=O)C2=CC=CC=C12",      # Anthraquinone
        "[O]=C1C=CC2=CC=CC=C2C1=O",          # Polycyclic-like
        "O=C1C=C2C3=CC=CC=C3C2=CC1=O",      # Tetracyclic-like
        "[O]=C1C2=C(C=CC=C2)C=CC1=O",       # Extended aromatic naphthalene
        "[O]=C1OC=CC2=C1C=CC=C2",           # Heterocyclic with oxygen
        "[C,N]1-2=[C-3]=[C,X2]=[C]=[C-2]=O" # Generalized quinoid ring with heteroatom inclusions
        # Include any other relevant quinone patterns as needed for broad coverage
    ]
    
    for pattern in quinone_patterns:
        try:
            smarts = Chem.MolFromSmarts(pattern)
            if smarts is None:
                continue
            if mol.HasSubstructMatch(smarts):
                return True, "Molecule contains a fully conjugated cyclic dione structure indicative of a quinone"
        except:
            continue

    return False, "Molecule does not contain a quinone structural feature"