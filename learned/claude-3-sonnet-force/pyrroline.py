"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:35447 pyrroline
Any organic heteromonocyclic compound with a structure based on a dihydropyrrole.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is an organic heteromonocyclic compound with a dihydropyrrole structure.

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
    
    # Check for pyrroline ring
    pyrroline_ring = Chem.MolFromSmarts("C1=NCCC1")
    pyrroline_matches = mol.GetSubstructMatches(pyrroline_ring)
    if not pyrroline_matches:
        return False, "No pyrroline ring found"
    
    # Enumerate tautomers and check for pyrroline structures
    tautomers = AllChem.EnumerateTautomers(mol)
    for tautomer in tautomers:
        tautomer_matches = tautomer.GetSubstructMatches(pyrroline_ring)
        if tautomer_matches:
            break
    else:
        return False, "No pyrroline tautomer found"
    
    # Check for disallowed heteroatoms and functional groups
    disallowed_smarts = ["[!#6&#7]", "*=O", "*=S", "*#*", "O=C=O"]
    disallowed_patterns = [Chem.MolFromSmarts(s) for s in disallowed_smarts]
    for pattern in disallowed_patterns:
        if mol.HasSubstructMatch(pattern):
            return False, "Contains disallowed heteroatoms or functional groups"
    
    # Check for allowed substituents (alkyl, alkenyl, aryl, halogen, amine, ether)
    allowed_smarts = ["C", "c", "X", "N", "O"]
    allowed_patterns = [Chem.MolFromSmarts("[" + s + "]") for s in allowed_smarts]
    pyrroline_atoms = [atom.GetIdx() for match in tautomer_matches for atom in mol.GetAtomWithIdx(match)]
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in pyrroline_atoms:
            allowed = False
            for pattern in allowed_patterns:
                if atom.IsInRingSize(5) and mol.HasSubstructMatch(pattern, atom.GetIdx()):
                    allowed = True
                    break
            if not allowed:
                return False, "Contains disallowed substituents"
    
    return True, "Molecule contains a pyrroline ring system"