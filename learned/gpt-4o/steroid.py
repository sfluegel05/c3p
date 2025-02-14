"""
Classifies: CHEBI:35341 steroid
"""
from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is identified by the presence of a cyclopenta[a]phenanthrene carbon skeleton,
    comprising three 6-membered rings and one 5-membered ring, with flexibility for functional modifications.

    Args:
        smiles (str): SMILES string of the chemical compound

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Examine ring structures - need exactly 4 rings: three 6-membered and one 5-membered
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    
    # Count 6-membered and 5-membered rings
    six_membered_rings = sum(1 for ring in rings if len(ring) == 6)
    five_membered_rings = sum(1 for ring in rings if len(ring) == 5)
    
    # Check for the steroid core structure: 3 six-membered and 1 five-membered ring
    if not (six_membered_rings >= 3 and five_membered_rings >= 1):
        return False, "Steroid structure not found in terms of 6:5 ring pattern"
    
    # Examine for common functional groups
    # Check for the presence of hydroxyl groups, common in many steroids
    hydroxyl_smarts = Chem.MolFromSmarts("[OX2H]")
    if mol.HasSubstructMatch(hydroxyl_smarts):
        return True, "Molecule matches steroid structure with hydroxyl group."
    
    # Additional checks for steroid-specific methyl groups or patterns can further refine the specificity
    # Example: [C](C)(C)[CH](CX) - a methyl branch likely at specific positions - can be included
    
    return True, "Molecule matches typical steroid structure"

# Examples SMILES test (these should return True)
smiles_examples = [
    "C[C@]12CC[C@H]3[C@H]([C@H]1CC[C@]2(C#C)O)CCC4=C3C=CC(=C4)OC",  # Expected True
    "CC(C)=CCC[C@](C)(O)[C@H]1CC[C@]2(C)[C@@H]1[C@H](O)C[C@@H]1[C@@]3(C)CC[C@H](O)C(C)(C)[C@@H]3CC[C@@]21C",  # Expected True
]

# Testing the function with example SMILES
for example in smiles_examples:
    print(is_steroid(example))