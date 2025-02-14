"""
Classifies: CHEBI:48953 cyclohexenones
"""
from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is a six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for a 6 membered alicyclic ring with a carbonyl and a double bond in the ring.
    #  '[C;R1;X4,X3]=[O]' : Carbonyl group, part of ring
    #  '[C;R1;X4,X3]=[C;R1;X4,X3]' : Double bond part of ring
    # '[C;R1;X4,X3][C;R1;X4,X3][C;R1;X4,X3][C;R1;X4,X3][C;R1;X4,X3][C;R1;X4,X3]' : six membered ring
    # the [C;R1;X4,X3] means a carbon atom part of a ring and 3 or 4 bonds
    cyclohexenone_pattern = Chem.MolFromSmarts("[C;R1;X4,X3](=[O])[C;R1;X4,X3]~[C;R1;X4,X3]=[C;R1;X4,X3]~[C;R1;X4,X3]~[C;R1;X4,X3]")
    if mol.HasSubstructMatch(cyclohexenone_pattern):
        return True, "Molecule is a cyclohexenone"
    
    # check for a carbonyl and a double bond on the same atom in a 6 membered ring
    # a double bond to the ring is represented by ~ 
    cyclohexenone_pattern_2 = Chem.MolFromSmarts("[C;R1;X3](=[O])=[C;R1;X4,X3]~[C;R1;X4,X3]~[C;R1;X4,X3]~[C;R1;X4,X3]~[C;R1;X4,X3]")
    if mol.HasSubstructMatch(cyclohexenone_pattern_2):
        return True, "Molecule is a cyclohexenone"
    
    return False, "Molecule does not contain a six-membered alicyclic ring with a carbonyl and double bond"