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

    #Check for a minimum number of heavy atoms (C, N, O, S, P etc)
    heavy_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1)
    if heavy_atom_count < 8:
       return False, "Too small, likely not alkaloid"

    # Check if there is at least one *heterocyclic* nitrogen - nitrogen that is part of the ring system
    heterocyclic_nitrogen_pattern = Chem.MolFromSmarts("[#7;R1;!H0]") # Must be part of a ring, and not attached to H
    has_heterocyclic_nitrogen = mol.HasSubstructMatch(heterocyclic_nitrogen_pattern)
    
    # Check for exocyclic nitrogen, excluding amides
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and not atom.IsInRing():
            
            is_amide_or_similar = False # Check for amides
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6:
                    for neighbor2 in neighbor.GetNeighbors():
                        if neighbor2.GetAtomicNum() == 8 and mol.GetBondBetweenAtoms(neighbor.GetIdx(), neighbor2.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                            is_amide_or_similar = True
                            break
                    if is_amide_or_similar:
                        break # Move to next N
            if not is_amide_or_similar:
                #Check if nitrogen is only attached to C atoms and no other non H
                all_c_neighbors = True
                carbon_count = 0
                for neighbor in atom.GetNeighbors():
                      if neighbor.GetAtomicNum() !=6 and neighbor.GetAtomicNum() != 1:
                         all_c_neighbors = False
                      if neighbor.GetAtomicNum() == 6:
                        carbon_count+=1

                if all_c_neighbors and carbon_count < 3:
                   return False, "Contains simple exocyclic amine"

    reason = "Contains at least one nitrogen atom, "
    if has_ring:
        reason += "has at least one ring, "
    if has_heterocyclic_nitrogen:
      reason += "and at least one N in a ring, "

    reason += "and is large enough, excluding simple exocyclic amine"


    return True, reason