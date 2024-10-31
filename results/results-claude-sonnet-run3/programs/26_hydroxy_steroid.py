from rdkit import Chem
from rdkit.Chem import AllChem

def is_26_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 26-hydroxy steroid.
    A 26-hydroxy steroid is any steroid with a hydroxy group at position 26.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 26-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Basic steroid core pattern (four fused rings)
    steroid_core = Chem.MolFromSmarts("[#6]1~[#6]~[#6]~[#6]2~[#6]~[#6]~[#6]3~[#6]~[#6]~[#6]4~[#6]~[#6]~[#6]~[#6]~[#6]4~[#6]~[#6]3~[#6]~[#6]2~[#6]~1")
    
    # Pattern for terminal CH2OH group with proper chain length for position 26
    # Match any terminal CH2OH group
    terminal_ch2oh = Chem.MolFromSmarts("[CH2][OH]")

    if not mol.HasSubstructMatch(steroid_core):
        return False, "Does not contain steroid core"

    if not mol.HasSubstructMatch(terminal_ch2oh):
        return False, "Does not contain terminal CH2OH group"

    # Count carbons and oxygens
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    oxygen_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8])
    
    if carbon_count < 21:  # Minimum carbons for a steroid
        return False, "Insufficient carbon count for steroid structure"

    # Get all terminal CH2OH groups
    ch2oh_matches = mol.GetSubstructMatches(terminal_ch2oh)
    
    # For each terminal CH2OH group, check if it's connected to the steroid core
    # through an appropriate length chain
    for match in ch2oh_matches:
        ch2_atom = mol.GetAtomWithIdx(match[0])
        # Get all atoms within 6 bonds of the CH2 group
        atoms_in_range = set()
        for atom in mol.GetAtoms():
            if atom.GetIdx() == ch2_atom.GetIdx():
                continue
            path = Chem.GetShortestPath(mol, ch2_atom.GetIdx(), atom.GetIdx())
            if path and len(path) <= 7:  # CH2OH should be about 6-7 bonds away from steroid core
                atoms_in_range.add(atom.GetIdx())
        
        # If we find a path of appropriate length to the steroid core
        if atoms_in_range:
            return True, "Contains steroid core and 26-hydroxy group"

    return False, "CH2OH group not properly positioned for 26-hydroxylation"
# Pr=None
# Recall=0.0