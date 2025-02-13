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

    # More flexible porphyrin pattern - just looking for 4 nitrogens connected to rings
    porphyrin_pattern = Chem.MolFromSmarts("[#7]~1~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#6]~[#7]~[#6]~[#6]~[#6]~[#6]~[#7]~1")
    if porphyrin_pattern is None:
        return None, "Invalid SMARTS pattern for porphyrin"
    if not mol.HasSubstructMatch(porphyrin_pattern):
        return False, "No porphyrin core structure found"

    # Look for fifth ring - more flexible pattern
    # Could be cyclopentanone or similar 5-membered ring
    fifth_ring_pattern = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]~[#6]1")
    if fifth_ring_pattern is None:
        return None, "Invalid SMARTS pattern for fifth ring"
    
    # Count number of 5-membered rings
    five_rings = len(mol.GetSubstructMatches(fifth_ring_pattern))
    if five_rings < 1:
        return False, "No fifth ring found"

    # Check for characteristic substituents
    substituents = []
    
    # Vinyl group
    vinyl_pattern = Chem.MolFromSmarts("C=C")
    if vinyl_pattern and mol.HasSubstructMatch(vinyl_pattern):
        substituents.append("vinyl")
        
    # Ethyl group
    ethyl_pattern = Chem.MolFromSmarts("CC")
    if ethyl_pattern and mol.HasSubstructMatch(ethyl_pattern):
        substituents.append("ethyl")
        
    # Methyl group
    methyl_pattern = Chem.MolFromSmarts("[CH3]")
    if methyl_pattern and mol.HasSubstructMatch(methyl_pattern):
        substituents.append("methyl")
        
    # Carboxyl/ester groups
    carboxyl_pattern = Chem.MolFromSmarts("[#6](=[#8])[#8]")
    if carboxyl_pattern and mol.HasSubstructMatch(carboxyl_pattern):
        substituents.append("carboxyl/ester")

    if not substituents:
        return False, "No characteristic substituents found"

    # Count carbons and nitrogens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if c_count < 20:
        return False, "Too few carbons for chlorophyll"
    if n_count != 4:
        return False, "Must have exactly 4 nitrogens in porphyrin core"

    # Additional check for conjugated system
    conjugated_pattern = Chem.MolFromSmarts("c:c:c:c:c")
    if conjugated_pattern and not mol.HasSubstructMatch(conjugated_pattern):
        return False, "Missing conjugated ring system"

    return True, f"Contains magnesium porphyrin core with fifth ring and {', '.join(substituents)} groups"