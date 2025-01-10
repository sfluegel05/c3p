"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: chlorophyll molecules
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    Chlorophylls are magnesium porphyrins with a fifth ring and various side chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chlorophyll, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for magnesium atom
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == 'Mg']
    if not mg_atoms:
        return False, "No magnesium atom found"
    if len(mg_atoms) > 1:
        return False, "More than one magnesium atom found"

    # Check for porphyrin-like core (4 nitrogens around Mg)
    mg = mg_atoms[0]
    neighboring_n = len([atom for atom in mg.GetNeighbors() if atom.GetSymbol() == 'N'])
    if neighboring_n != 4:
        return False, f"Magnesium should have 4 coordinating nitrogens, found {neighboring_n}"

    # Look for the basic porphyrin structure with 4 pyrrole-like rings
    porphyrin_pattern = Chem.MolFromSmarts("[#7]1:[#6]:[#6]:[#6]2:[#6]:[#6]:[#6]:[#6]:[#7]:[#6]:[#6]:[#6]3:[#6]:[#6]:[#6]:[#6]:[#7]:[#6]:[#6]:[#6]4:[#6]:[#6]:[#6]:[#6]:[#7]:1:[Mg]2134")
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin core structure found"

    # Look for fifth ring (cyclopentanone ring)
    cyclopentanone_pattern = Chem.MolFromSmarts("C1CC(=O)C(C1)=O")
    if not mol.HasSubstructMatch(cyclopentanone_pattern):
        return False, "No fifth ring (cyclopentanone) found"

    # Check for common substituents found in chlorophylls
    # Vinyl group
    vinyl_pattern = Chem.MolFromSmarts("C=C")
    # Ethyl group
    ethyl_pattern = Chem.MolFromSmarts("CC")
    # Methyl group
    methyl_pattern = Chem.MolFromSmarts("C[CH3]")
    
    substituents_found = []
    if mol.HasSubstructMatch(vinyl_pattern):
        substituents_found.append("vinyl")
    if mol.HasSubstructMatch(ethyl_pattern):
        substituents_found.append("ethyl")
    if mol.HasSubstructMatch(methyl_pattern):
        substituents_found.append("methyl")
        
    if not substituents_found:
        return False, "No characteristic substituents found"

    # Look for ester groups (common in chlorophylls)
    ester_pattern = Chem.MolFromSmarts("[#6](=O)-O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester groups found"

    # Count carbons and nitrogens to ensure reasonable size
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20:
        return False, "Too few carbons for chlorophyll"
    if n_count != 4:
        return False, "Must have exactly 4 nitrogens in porphyrin core"

    return True, f"Contains magnesium porphyrin core with fifth ring and {', '.join(substituents_found)} substituents"