"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Basic steroid skeleton patterns (multiple SMARTS to account for different oxidation states)
    steroid_patterns = [
        # Basic steroid skeleton with saturated rings
        "[#6]12[#6][#6][#6]3[#6][#6][#6]4[#6][#6][#6][#6]([#6]1)[#6][#6][#6]2[#6]34",
        # Steroid skeleton with possible double bonds
        "[#6]12[#6][#6][#6]3[#6][#6][#6]4[#6]=,:[#6][#6]=,:[#6]([#6]1)[#6][#6][#6]2[#6]34",
    ]
    
    found_skeleton = False
    for pattern in steroid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_skeleton = True
            break
            
    if not found_skeleton:
        return False, "No steroid skeleton found"

    # Count rings
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 4:
        return False, "Insufficient number of rings for steroid structure"

    # Check carbon count (steroids typically have 17+ carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 17:
        return False, "Too few carbons for steroid structure"

    # Check for typical molecular weight range of steroids (240-650 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if not (240 <= mol_wt <= 650):
        return False, f"Molecular weight {mol_wt:.1f} outside typical steroid range"

    # Check for typical steroid characteristics
    # Angular methyl groups pattern (at C-10 and C-13)
    methyl_pattern = Chem.MolFromSmarts("[CH3]C12[CH2,CH][CH2,CH][#6]3[#6][#6][#6]4[#6][#6][#6][#6]([#6]1)[#6][#6][#6]2[#6]34")
    if not mol.HasSubstructMatch(methyl_pattern):
        # Some steroids might have lost methyl groups, so this is not a strict requirement
        pass

    # Look for typical steroid functional groups
    functional_groups = {
        "hydroxyl": "[OH]",
        "ketone": "C(=O)C",
        "ester": "C(=O)O[#6]",
        "carboxylic acid": "C(=O)O[H]",
    }
    
    found_groups = []
    for group_name, smarts in functional_groups.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            found_groups.append(group_name)

    # Final classification
    if found_skeleton:
        reason = "Contains steroid skeleton"
        if found_groups:
            reason += f" with {', '.join(found_groups)} groups"
        return True, reason
    
    return False, "Structure does not match steroid characteristics"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:35341',
        'name': 'steroid',
        'definition': 'Any of naturally occurring compounds and synthetic analogues, '
                     'based on the cyclopenta[a]phenanthrene carbon skeleton.',
        'parents': ['CHEBI:33860']
    }
}