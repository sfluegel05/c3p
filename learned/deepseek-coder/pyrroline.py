"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:38943 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is a dihydropyrrole ring, which is a five-membered ring with one nitrogen and one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more general pyrroline pattern: a five-membered ring with one nitrogen and one double bond
    pyrroline_pattern = Chem.MolFromSmarts("[nH0]1[CH]=,:[CH][CH2]1")
    pyrroline_pattern_alt = Chem.MolFromSmarts("[nH0]1[CH2][CH]=,:[CH]1")
    pyrroline_pattern_alt2 = Chem.MolFromSmarts("[nH0]1[CH]=,:[CH2][CH]1")
    
    # Check if the molecule contains any of the pyrroline patterns
    if (mol.HasSubstructMatch(pyrroline_pattern) or 
        mol.HasSubstructMatch(pyrroline_pattern_alt) or 
        mol.HasSubstructMatch(pyrroline_pattern_alt2)):
        
        # Ensure the matched ring is a five-membered ring
        ring_info = mol.GetRingInfo()
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 7:  # Nitrogen atom
                for ring in ring_info.AtomRings():
                    if len(ring) == 5 and atom.GetIdx() in ring:
                        return True, "Contains a dihydropyrrole ring (pyrroline)"
        
    return False, "No dihydropyrrole ring (pyrroline) found"