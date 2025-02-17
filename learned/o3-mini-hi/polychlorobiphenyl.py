"""
Classifies: CHEBI:53156 polychlorobiphenyl
"""
"""
Classifies: Polychlorobiphenyl
Definition: A biphenyl compound containing between 2 and 10 chlorine atoms attached to the two benzene rings.
Additional Requirement: The only substituents directly attached to atoms in the biphenyl (i.e. the two benzene rings) must be chlorine atoms or hydrogens.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_polychlorobiphenyl(smiles: str):
    """
    Determines if a molecule is a polychlorobiphenyl based on its SMILES string.
    
    Requirements:
      - The molecule must contain a biphenyl scaffold: two benzene rings (six-membered aromatic rings composed solely of carbons)
        connected by a single bond.
      - For the atoms that belong to these two rings, the only substituents allowed (aside from hydrogens) are chlorine atoms.
      - The total number of chlorine atoms attached to the core rings must be between 2 and 10.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings in the molecule.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Identify benzene rings: 6-membered, all carbons, and aromatic.
    benzene_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            # Check that each atom in ring is aromatic carbon
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and mol.GetAtomWithIdx(idx).GetIsAromatic() 
                   for idx in ring):
                benzene_rings.append(set(ring))
    
    if len(benzene_rings) < 2:
        return False, "Less than two benzene rings found"
    
    # Look for a pair of benzene rings connected by a single bond.
    biphenyl_pairs = []
    for i in range(len(benzene_rings)):
        for j in range(i+1, len(benzene_rings)):
            found_connection = False
            # Check if any atom from ring i is connected to any atom from ring j
            for idx in benzene_rings[i]:
                atom = mol.GetAtomWithIdx(idx)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in benzene_rings[j]:
                        bond = mol.GetBondBetweenAtoms(idx, neighbor.GetIdx())
                        # We require a single bond connecting the two rings
                        if bond is not None and bond.GetBondType() == rdchem.BondType.SINGLE:
                            found_connection = True
                            break
                if found_connection:
                    break
            if found_connection:
                biphenyl_pairs.append((benzene_rings[i], benzene_rings[j]))
    
    if not biphenyl_pairs:
        return False, "No biphenyl scaffold (two connected benzene rings) found"
    
    # For each candidate biphenyl pair, check the substituents on the ring atoms.
    for ring1, ring2 in biphenyl_pairs:
        biphenyl_atoms = ring1.union(ring2)
        chlorine_count = 0
        valid_core = True
        disallowed = None
        
        # For each atom in the biphenyl, inspect its neighbors that are not part of the ring.
        for idx in biphenyl_atoms:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in biphenyl_atoms:
                    continue  # skip atoms in the core
                # Allowed: hydrogen (atomic num 1) or chlorine (atomic num 17)
                anum = neighbor.GetAtomicNum()
                if anum == 17:
                    chlorine_count += 1
                elif anum == 1:
                    continue
                else:
                    valid_core = False
                    disallowed = neighbor.GetSymbol()
                    # Once a disallowed substituent is found, we break out
                    break
            if not valid_core:
                break
        
        # Now check if this biphenyl pair qualifies.
        if not valid_core:
            # For this candidate, we encountered a disallowed substituent.
            # We do not immediately return False in order to try other biphenyl pairs.
            continue
        
        # Check the chlorine count range.
        if chlorine_count < 2:
            continue
        if chlorine_count > 10:
            continue
        
        # If we get here, we have a valid biphenyl satisfying all conditions.
        return True, f"Contains a biphenyl scaffold with {chlorine_count} chlorine substituents and no extra groups on the rings"
    
    # If none of the candidate biphenyl pairs qualified, report the issue.
    return False, "No biphenyl scaffold with the required chlorine substitution pattern found"

# Example calls (for testing purposes):
if __name__ == "__main__":
    test_examples = [
        ("Clc1cc(Cl)c(Cl)c(c1)-c1cc(Cl)cc(Cl)c1Cl", "2,2',3,3',5,5'-hexachlorobiphenyl (should be accepted)"),
        ("Clc1ccc(Cl)c(c1)-c1ccccc1", "2,5-dichlorobiphenyl (should be accepted)"),
        ("Oc1c(Cl)cc(Cl)cc1-c1ccccc1", "A biphenyl with an -OH (should be rejected)")
    ]
    for smi, desc in test_examples:
        result, reason = is_polychlorobiphenyl(smi)
        print(f"{desc}\n  SMILES: {smi}\n  Result: {result}, Reason: {reason}\n")