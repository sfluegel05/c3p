"""
Classifies: CHEBI:35381 monosaccharide
"""
from rdkit import Chem

def is_monosaccharide(smiles: str):
    """
    Determines if a molecule is a monosaccharide based on its SMILES string.
    Monosaccharides are polyhydroxy aldehydes or ketones, with at least three carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monosaccharide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")
    
    # Identify the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 3:
        return (False, "Too few carbon atoms; monosaccharides must have at least 3")
    
    # Identify the number of hydroxyl groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    # Hydroxyls would be counted based on whether they are -OH and not involved in saccharide links
    hydroxyl_count = len([1 for bond in mol.GetBonds() 
                          if bond.GetBeginAtom().GetAtomicNum() == 8 or bond.GetEndAtom().GetAtomicNum() == 8])
    if hydroxyl_count < 2:
        return (False, f"Insufficient hydroxyl groups found; needed at least 2, found {hydroxyl_count}")

    # Check for potential carbonyl group, which might be oxidized or involved in ring formation
    # Even for cyclic forms assume possible intramolecular hemiacetal for saccharide potential
    possible_hemiacetal = False
    for atom in mol.GetAtoms():
        if (atom.GetAtomicNum() == 6 and 
            any(neigh.GetAtomicNum() == 8 for neigh in atom.GetNeighbors())):
            if (len([neigh for neigh in atom.GetNeighbors() if neigh.GetAtomicNum() == 1]) > 0):
                possible_hemiacetal = True
                
    if not possible_hemiacetal:
        return (False, "No potential carbonyl functionality detected.")
    
    # Verify molecular weight and ensure it does not exceed typical monosaccharide weight (~200 Da)
    mol_weight = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight > 300:
        return (False, "Molecular weight too high; suggests multi-unit composition.")

    return (True, "Structure is a monosaccharide with cyclic/acyclic carbonyl potential and sufficient hydroxyls.")