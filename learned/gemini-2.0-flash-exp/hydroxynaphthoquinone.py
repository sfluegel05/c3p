"""
Classifies: CHEBI:132155 hydroxynaphthoquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hydroxynaphthoquinone(smiles: str):
    """
    Determines if a molecule is a hydroxynaphthoquinone based on its SMILES string.
    A hydroxynaphthoquinone is a naphthoquinone with at least one hydroxy substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxynaphthoquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for naphthoquinone core and hydroxy
    naphthoquinone_pattern1 = Chem.MolFromSmarts("C1=CC=C2C(=O)C=CC(=O)C2=C1")
    naphthoquinone_pattern2 = Chem.MolFromSmarts("C1=CC=CC2=C1C(=O)C=CC2=O")
    hydroxy_pattern = Chem.MolFromSmarts("[OH]")
    
    # Check if core is present
    if not (mol.HasSubstructMatch(naphthoquinone_pattern1) or mol.HasSubstructMatch(naphthoquinone_pattern2)):
        return False, "No naphthoquinone core found"
    
    # Check if at least one hydroxy is directly attached to core
    matches = mol.GetSubstructMatches(hydroxy_pattern)

    hydroxylated = False
    for match in matches:
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if (mol.HasSubstructMatch(naphthoquinone_pattern1, useAllAtoms=True) and neighbor.GetIdx() in mol.GetSubstructMatch(naphthoquinone_pattern1)) or \
                    (mol.HasSubstructMatch(naphthoquinone_pattern2, useAllAtoms=True) and neighbor.GetIdx() in mol.GetSubstructMatch(naphthoquinone_pattern2)) :
                        hydroxylated = True
                        break
            if hydroxylated:
                break
        if hydroxylated:
            break
            
    if not hydroxylated:
         return False, "No hydroxy group directly attached to naphthoquinone core"

    return True, "Contains a naphthoquinone core with at least one hydroxy substituent."