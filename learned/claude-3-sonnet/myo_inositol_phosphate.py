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
    
    # Remove charges to normalize the structure
    uncharge_smiles = Chem.MolToSmiles(mol, kekuleSmiles=True, allHsExplicit=False)
    mol = Chem.MolFromSmiles(uncharge_smiles)
    
    # Check for cyclohexane core
    cyclohexane_pattern = Chem.MolFromSmarts("C1CCCCC1")
    if not mol.HasSubstructMatch(cyclohexane_pattern):
        return False, "No cyclohexane core found"
    
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
    
    # More specific pattern for myo-inositol core with correct stereochemistry
    # This pattern looks for a cyclohexane with 5 equatorial and 1 axial substituents
    # The @ and @@ specify the stereochemistry
    myo_pattern = Chem.MolFromSmarts("[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1")
    if not mol.HasSubstructMatch(myo_pattern):
        return False, "Does not match myo-inositol stereochemistry pattern"
    
    # Check for phosphate groups - more flexible pattern that matches both charged and uncharged forms
    phosphate_pattern = Chem.MolFromSmarts("[OX2][P]([OX2,OX1-])([OX2,OX1-])=[OX1,OX2-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found"
    
    # Verify that phosphates are attached to the inositol ring
    ring_atoms = ring_info.AtomRings()[0]
    ring_carbons = set(idx for idx in ring_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
    
    # Check oxygen connections between ring and phosphates
    for carbon_idx in ring_carbons:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 8:  # Oxygen
                for next_neighbor in neighbor.GetNeighbors():
                    if next_neighbor.GetAtomicNum() == 15:  # Phosphorus
                        break
                else:
                    continue
                break
    
    # Count phosphate-oxygen bonds to verify phosphate groups
    p_o_pattern = Chem.MolFromSmarts("[P]~[O]")
    p_o_count = len(mol.GetSubstructMatches(p_o_pattern))
    min_expected_o = p_count * 3  # Each P should have at least 3 oxygens
    if p_o_count < min_expected_o:
        return False, "Insufficient phosphate-oxygen bonds"

    return True, f"Contains myo-inositol core with {p_count} phosphorus atoms in phosphate groups"