"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for cyclohexane core (more general pattern)
    cyclohexane_pattern = Chem.MolFromSmarts("C1CCCCC1")
    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "No cyclohexane core found"
    
    # Check for phosphate groups (including mono-, di-, and tri-phosphates)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P](=[OX1])([OX2,OX1-])[OX2,OX1-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found"
    
    # Count carbons (should be exactly 6 for myo-inositol)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 6:
        return False, f"Found {c_count} carbons, need exactly 6 for myo-inositol"
    
    # Check for phosphorus
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count == 0:
        return False, "No phosphorus atoms found"
        
    # Check for ring size (should be 6)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    if not any(size == 6 for size in ring_sizes):
        return False, "No 6-membered ring found"
        
    # Check that all carbons in the ring have oxygen substituents
    # This includes both OH and OP groups
    ring_atoms = ring_info.AtomRings()[0]  # Get atoms in the first 6-membered ring
    for ring_atom_idx in ring_atoms:
        atom = mol.GetAtomWithIdx(ring_atom_idx)
        if atom.GetAtomicNum() != 6:  # Skip if not carbon
            continue
        
        # Check if carbon has oxygen neighbor
        has_oxygen = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                has_oxygen = True
                break
        if not has_oxygen:
            return False, "Not all ring carbons have oxygen substituents"
    
    # Check for characteristic myo-inositol substitution pattern
    # Look for 6-membered ring with all carbons having one oxygen substituent
    myo_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C(O)C(O)1")
    if not mol.HasSubstructMatch(myo_pattern):
        return False, "Does not match myo-inositol substitution pattern"
    
    # Count total phosphate-related oxygens
    # Each P should have at least 3 oxygens (mono-phosphate)
    p_o_pattern = Chem.MolFromSmarts("[P]~[O]")
    p_o_count = len(mol.GetSubstructMatches(p_o_pattern))
    if p_o_count < p_count * 3:
        return False, "Insufficient phosphate-oxygen bonds"

    return True, f"Contains myo-inositol core with {p_count} phosphorus atoms in phosphate groups"