"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: Amino sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar with one or more hydroxyl groups replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for 5 and 6 membered rings with one oxygen and carbons
    sugar_ring_patterns = [
        Chem.MolFromSmarts("[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3]1"), #5 membered ring
        Chem.MolFromSmarts("[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3][CX4,CX3]1")  # 6 membered ring
    ]

    has_sugar_ring = False
    ring_atoms = []
    for pattern in sugar_ring_patterns:
        if mol.HasSubstructMatch(pattern):
          has_sugar_ring = True
          matches = mol.GetSubstructMatches(pattern)
          for match in matches:
                ring_atoms.extend(match)
          break #take only the first ring

    if not has_sugar_ring:
        return False, "No sugar ring structure found"


    #check the atoms of the ring. there must be at most 1 oxygen, all other carbon and all carbons must be saturated
    oxygen_count = 0
    carbon_count = 0
    for atom_index in ring_atoms:
          atom = mol.GetAtomWithIdx(atom_index)
          if atom.GetAtomicNum() == 8:
              oxygen_count += 1
          elif atom.GetAtomicNum() == 6 and atom.GetHybridization() != Chem.HybridizationType.SP2: #exclude not saturated carbons
                carbon_count += 1
          else:
               return False, "Ring must contain only carbons and one oxygen"


    if oxygen_count != 1:
        return False, "Sugar ring must have exactly one oxygen"

    # Check for amino groups (NH2, NHR, NR2) and also acetamido groups attached directly or indirectly to a ring carbon
    amino_group_pattern = Chem.MolFromSmarts("[NX3;!H0]-[CX4,CX3]") #amino directly bonded to a C
    substituted_amino_pattern = Chem.MolFromSmarts("[NX3;!H0][CX3](=[OX1])") #acetamide -NH-C=O
    amino_group_pattern_indirect = Chem.MolFromSmarts("[NX3;!H0]~[CX4,CX3]~[CX4,CX3]1[OX2][CX4,CX3][CX4,CX3][CX4,CX3]1") #amino group via linker to ring
    
    
    amino_matches = []
    for atom_index in ring_atoms:
        atom = mol.GetAtomWithIdx(atom_index)
        if atom.GetAtomicNum() == 6: #check only for carbons of the ring
            for neighbor in atom.GetNeighbors():
                  if neighbor.Match(amino_group_pattern) or neighbor.Match(substituted_amino_pattern):
                      amino_matches.append(neighbor)
                      
    amino_matches2 = mol.GetSubstructMatches(amino_group_pattern_indirect)
    if not amino_matches and not amino_matches2:
       return False, "No amino group substitution found"

    # If both conditions are met, it is an amino sugar
    return True, "Contains a sugar ring with one or more amino group substitutions."