"""
Classifies: CHEBI:22315 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is likely an alkaloid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely an alkaloid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one nitrogen atom
    nitrogen_pattern = Chem.MolFromSmarts("[#7]")
    if not mol.HasSubstructMatch(nitrogen_pattern):
        return False, "No nitrogen atom present"
    
    # Check for presence of at least one ring
    ring_pattern = Chem.MolFromSmarts("[R]")
    has_ring = mol.HasSubstructMatch(ring_pattern)

    # Check if there is at least one nitrogen in a ring - this is not a strict condition
    heterocyclic_nitrogen_pattern = Chem.MolFromSmarts("[#7;R]")
    has_ring_nitrogen = mol.HasSubstructMatch(heterocyclic_nitrogen_pattern)
    
    #Check for a minimum number of C and heteroatoms (N, O, S)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    hetero_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in [7, 8, 16])

    if c_count < 10 or hetero_count < 2:
        return False, "Too small, likely not alkaloid"


    # Check for exocyclic nitrogen, excluding amides
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and not atom.IsInRing():
            
            is_amide_or_similar = False
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    for neighbor2 in neighbor.GetNeighbors():
                        if neighbor2.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(neighbor.GetIdx(), neighbor2.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                            is_amide_or_similar = True
                            break
                    if is_amide_or_similar:
                        break # Move to next N
            if not is_amide_or_similar:
                #Check if nitrogen is only attached to C atoms. If yes, it is considered an exocyclic amine
                neighbor_carbon_count = 0
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6:
                        neighbor_carbon_count +=1
                all_c_neighbors = True
                for neighbor in atom.GetNeighbors():
                      if neighbor.GetAtomicNum() !=6:
                         all_c_neighbors = False

                if all_c_neighbors and neighbor_carbon_count < 3: # simple amines should not have 3 C neighbors
                   return False, "Contains simple exocyclic amine"



    reason = "Contains at least one nitrogen atom, "
    if has_ring:
        reason += "has at least one ring, "
    if has_ring_nitrogen:
      reason += "and at least one N in a ring, "

    reason += "and is large enough, excluding simple exocyclic amine"


    return True, reason