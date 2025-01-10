"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: CHEBI:36275 lactol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    A lactol is a cyclic hemiacetal formed by intramolecular addition of a hydroxy group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Verify that the molecule is cyclic
    if not mol.GetRingInfo().NumRings():
        return False, "Molecule is not cyclic (lactols are cyclic hemiacetals)"

    # Look for the general lactol pattern: a ring with an oxygen and a hydroxyl group on the same carbon
    lactol_pattern = Chem.MolFromSmarts("[C;R][OH][C;R][O;R]")
    if not mol.HasSubstructMatch(lactol_pattern):
        return False, "No lactol pattern found (cyclic hemiacetal structure)"

    # Check for the presence of a carbonyl group (aldehyde or ketone) in the molecule
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(carbonyl_pattern):
        return False, "No carbonyl group found (aldehyde or ketone required for lactol formation)"

    # Check that the hydroxyl group and oxygen are in the same ring
    ring_info = mol.GetRingInfo()
    for match in mol.GetSubstructMatches(lactol_pattern):
        carbon1 = match[0]
        hydroxyl = match[1]
        carbon2 = match[2]
        oxygen = match[3]
        
        # Verify all atoms are in the same ring
        if not any(ring for ring in ring_info.AtomRings() 
                  if carbon1 in ring and hydroxyl in ring and 
                     carbon2 in ring and oxygen in ring):
            return False, "Hydroxyl group and oxygen not in the same ring system"

        # Verify the hydroxyl is on the same carbon as the ring oxygen
        if mol.GetBondBetweenAtoms(carbon2, oxygen) is None:
            return False, "Hydroxyl group not on same carbon as ring oxygen"

    # Check that the carbonyl group is in the same ring system
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    for carbonyl_match in carbonyl_matches:
        carbonyl_atom = carbonyl_match[0]
        for ring in ring_info.AtomRings():
            if carbonyl_atom in ring:
                # Check if the carbonyl is in the same ring system as the lactol
                for match in mol.GetSubstructMatches(lactol_pattern):
                    if any(atom in ring for atom in match):
                        return True, "Contains a cyclic hemiacetal structure with a hydroxyl group and an oxygen in the ring"

    return False, "No carbonyl group in the same ring system as the lactol structure"

# Example usage:
# smiles = "C1C(C(C(O1)O)O)O"  # 2,3,4-trihydroxytetrahydrofuran
# print(is_lactol(smiles))