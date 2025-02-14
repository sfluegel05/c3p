"""
Classifies: CHEBI:65323 monoterpenoid indole alkaloid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_monoterpenoid_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is a monoterpenoid indole alkaloid based on its SMILES string.
    A monoterpenoid indole alkaloid contains an indole moiety, a monoterpenoid unit and is an alkaloid.
    It is biosynthesised from L-tryptophan and diisoprenoid (usually secolaganin) building blocks.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoterpenoid indole alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for a substituted Indole/Tryptamine-like substructure fused to a ring
    #   This pattern includes a substituted indole connected to a 5 or 6-membered ring
    #   The connection can be a direct bond or a spiro connection.
    #   The nitrogen can be substituted or charged.
    indole_pattern1 = Chem.MolFromSmarts("[c]1[c]([#6])[n][c][c]2[c]1[#6][#6]2") # Neutral N
    indole_pattern2 = Chem.MolFromSmarts("[c]1[c]([#6])[N][c][c]2[c]1[#6][#6]2") # Neutral N
    indole_pattern3 = Chem.MolFromSmarts("[c]1[c]([#6])[n+][c][c]2[c]1[#6][#6]2") # Charged N
    indole_pattern4 = Chem.MolFromSmarts("[c]1[c]([#6])[nH+][c][c]2[c]1[#6][#6]2") # Protonated N

    has_indole = mol.HasSubstructMatch(indole_pattern1) or mol.HasSubstructMatch(indole_pattern2) or mol.HasSubstructMatch(indole_pattern3) or mol.HasSubstructMatch(indole_pattern4)
    if not has_indole:
         return False, "No appropriate substituted indole substructure fused to a ring found."
    
    # 2. Check for a monoterpenoid-like substructure
    # Monoterpenoids are C10 units and often include oxygen bridges or cyclic systems or a di-isoprenoid structure.
    # Here we can try a few different ones for better coverage

    monoterpenoid_pattern1 = Chem.MolFromSmarts("[CX4]1[CX4][OX2][CX4]([CX4])([CX4])1") # Oxygen bridge in a ring
    monoterpenoid_pattern2 = Chem.MolFromSmarts("[CX4]1[CX4][CX4][CX4]([CX4])([CX4])1") # Cyclic pattern
    monoterpenoid_pattern3 = Chem.MolFromSmarts("[CX4]1([CX4])([CX4])[CX4][CX4]([CX4])([CX4])1") # another cyclic pattern
    monoterpenoid_pattern4 = Chem.MolFromSmarts("[CX4]1[CX4][CX4]=[CX3][CX4]([CX4])1") # unsaturation
    
    has_monoterpenoid = mol.HasSubstructMatch(monoterpenoid_pattern1) or mol.HasSubstructMatch(monoterpenoid_pattern2) or mol.HasSubstructMatch(monoterpenoid_pattern3) or mol.HasSubstructMatch(monoterpenoid_pattern4)
    if not has_monoterpenoid:
        return False, "No monoterpenoid-like substructure found"


    # 3. Check for Connection: The indole must be connected to the monoterpenoid
    # Simplified connectivity search.
    connected_pattern1 = Chem.MolFromSmarts("([c]1[c]([#6])[n][c][c]2[c]1[#6][#6]2)~[#6]~[CX4]1[CX4][OX2][CX4]([CX4])([CX4])1")
    connected_pattern2 = Chem.MolFromSmarts("([c]1[c]([#6])[n][c][c]2[c]1[#6][#6]2)~[#6]~[CX4]1[CX4][CX4][CX4]([CX4])([CX4])1")
    connected_pattern3 = Chem.MolFromSmarts("([c]1[c]([#6])[N][c][c]2[c]1[#6][#6]2)~[#6]~[CX4]1[CX4][OX2][CX4]([CX4])([CX4])1")
    connected_pattern4 = Chem.MolFromSmarts("([c]1[c]([#6])[N][c][c]2[c]1[#6][#6]2)~[#6]~[CX4]1[CX4][CX4][CX4]([CX4])([CX4])1")
    connected_pattern5 = Chem.MolFromSmarts("([c]1[c]([#6])[n+][c][c]2[c]1[#6][#6]2)~[#6]~[CX4]1[CX4][OX2][CX4]([CX4])([CX4])1")
    connected_pattern6 = Chem.MolFromSmarts("([c]1[c]([#6])[n+][c][c]2[c]1[#6][#6]2)~[#6]~[CX4]1[CX4][CX4][CX4]([CX4])([CX4])1")
    connected_pattern7 = Chem.MolFromSmarts("([c]1[c]([#6])[nH+][c][c]2[c]1[#6][#6]2)~[#6]~[CX4]1[CX4][OX2][CX4]([CX4])([CX4])1")
    connected_pattern8 = Chem.MolFromSmarts("([c]1[c]([#6])[nH+][c][c]2[c]1[#6][#6]2)~[#6]~[CX4]1[CX4][CX4][CX4]([CX4])([CX4])1")
    
    if not (mol.HasSubstructMatch(connected_pattern1) or mol.HasSubstructMatch(connected_pattern2) or
            mol.HasSubstructMatch(connected_pattern3) or mol.HasSubstructMatch(connected_pattern4) or
            mol.HasSubstructMatch(connected_pattern5) or mol.HasSubstructMatch(connected_pattern6) or
            mol.HasSubstructMatch(connected_pattern7) or mol.HasSubstructMatch(connected_pattern8)):
             return False, "Indole and monoterpenoid substructures are not connected."

    # 4. Complexity Checks
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 3:
         return False, "Too few rings for a monoterpenoid indole alkaloid."


    # Count carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:
        return False, "Too few carbons for a monoterpenoid indole alkaloid."


    # 5. Check for an alkaloid nitrogen
    # check that at least one nitrogen is in a ring
    nitrogen_in_ring = Chem.MolFromSmarts("[N;R]") # nitrogen in a ring
    if not mol.HasSubstructMatch(nitrogen_in_ring):
        return False, "Must have at least one nitrogen in a ring (alkaloid nature)"

    # Check molecular weight (less restrictive)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
         return False, "Molecular weight is too low for a monoterpenoid indole alkaloid."

    return True, "Contains an appropriate indole substructure fused to a ring system and connected to a monoterpenoid-like substructure, and sufficient complexity"