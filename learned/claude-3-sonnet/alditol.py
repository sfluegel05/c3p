"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:15972 alditol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for rings in the core structure
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        # Check if the rings are part of the core structure or just substituents
        core_rings = False
        for atom in mol.GetAtoms():
            if atom.IsInRing() and atom.GetAtomicNum() == 6:
                # If carbon is in ring, it's likely a cyclic sugar
                core_rings = True
                break
        if core_rings:
            return False, "Core structure contains rings (not acyclic)"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 3:
        return False, "Too few carbons for an alditol (minimum 3)"
    if o_count < 3:
        return False, "Too few oxygens for an alditol (minimum 3)"

    # Look for continuous carbon chain with hydroxyls
    # Pattern for carbon chain with alternating OH groups
    alditol_chain_pattern = Chem.MolFromSmarts("[CH2X4][OX2H1].[CHX4]([OX2H1])[CHX4]([OX2H1])")
    if not mol.HasSubstructMatch(alditol_chain_pattern):
        return False, "Missing required alditol carbon chain pattern"

    # Check for terminal CH2OH groups
    terminal_oh_pattern = Chem.MolFromSmarts("[CH2X4][OX2H1]")
    terminal_oh_matches = mol.GetSubstructMatches(terminal_oh_pattern)
    
    if len(terminal_oh_matches) < 1:
        return False, "Missing terminal CH2OH group"

    # Calculate the longest carbon chain
    longest_chain = rdMolDescriptors.CalcMolFormula(mol).count('C')
    
    # For derivatives, check if modifications maintain alditol core
    non_cho_atoms = sum(1 for atom in mol.GetAtoms() 
                       if atom.GetAtomicNum() not in [1, 6, 8])
    
    if non_cho_atoms > 0:
        # Check if the modifications are acceptable for alditol derivatives
        # (e.g., acetyl groups, phosphates, etc.)
        acceptable_mods = ["[C](=O)[C]", "[P](=O)([O])[O]"]
        has_acceptable_mods = False
        for mod in acceptable_mods:
            pattern = Chem.MolFromSmarts(mod)
            if pattern and mol.HasSubstructMatch(pattern):
                has_acceptable_mods = True
                break
        
        if has_acceptable_mods:
            return True, "Valid alditol derivative with acceptable modifications"
        else:
            return False, "Contains non-standard modifications for alditol"

    # For pure alditols, verify hydroxyl-to-carbon ratio
    oh_groups = len(mol.GetSubstructMatches(Chem.MolFromSmarts("[OX2H1]")))
    oh_ratio = oh_groups / c_count
    
    if 0.8 <= oh_ratio <= 1.2 and longest_chain >= 3:
        return True, "Pure alditol structure detected"
    
    return False, "Does not match alditol structure requirements"