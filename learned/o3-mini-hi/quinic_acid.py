"""
Classifies: CHEBI:26493 quinic acid
"""
#!/usr/bin/env python
"""
Classifies: A cyclitol carboxylic acid (quinic acid and derivatives)

A quinic acid core is defined here as:
  - A six-membered (cyclohexane) ring composed solely of sp³-hybridized, nonaromatic carbons 
    with all bonds between ring atoms being single.
  - At least one of the ring carbons must be substituted with a free carboxyl group.
  
A free carboxyl group is determined by:
  - Being a carbon atom (atomic number 6) with exactly three neighbors.
  - Exactly two of these neighbors are oxygen atoms, with one attached via a double bond 
    (C=O) and the other attached via a single bond, where that oxygen bears at least one hydrogen 
    or a negative formal charge (indicating a free –COOH group rather than an ester).
  
This more relaxed approach avoids the pitfalls of our earlier strict substituent count.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_quinic_acid(smiles: str):
    """
    Determines if a molecule is a quinic acid derivative based on its SMILES string.
    
    The criteria for a quinic acid core are:
      1. The molecule must contain at least one six-membered ring (cyclohexane) in which every atom 
         is a nonaromatic sp3-hybridized carbon and every bond between adjacent ring atoms is single.
      2. At least one of the ring atoms must have a free carboxyl group attached directly.
         A free carboxyl group is identified as a carbon binder with exactly three neighbors 
         (the ring atom plus two oxygen atoms), where one oxygen forms a double bond and the other 
         a single bond with a hydrogen (or negative charge).
         
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a quinic acid derivative, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string and add explicit hydrogens to catch O-H groups.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Helper function to check if a given carbon atom qualifies as a free carboxyl group.
    def is_free_carboxyl(atom):
        # Must be carbon
        if atom.GetAtomicNum() != 6:
            return False
        # Expect exactly three neighbors: one from the ring and two oxygens.
        if atom.GetDegree() != 3:
            return False
        oxy_count = 0
        double_oxygen = False
        single_oxygen_with_H = False
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                oxy_count += 1
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    double_oxygen = True
                elif bond.GetBondType() == rdchem.BondType.SINGLE:
                    # Check if this oxygen has an explicit hydrogen or is deprotonated
                    for onbr in nbr.GetNeighbors():
                        if onbr.GetAtomicNum() == 1:
                            single_oxygen_with_H = True
                            break
                    if nbr.GetFormalCharge() < 0:
                        single_oxygen_with_H = True
        return (oxy_count == 2 and double_oxygen and single_oxygen_with_H)
    
    # Get the ring information from the molecule.
    rings = mol.GetRingInfo().AtomRings()
    
    # Loop through candidate rings.
    for ring in rings:
        # Consider only six-membered rings.
        if len(ring) != 6:
            continue
            
        valid_ring = True
        # Check that every atom in the ring is a nonaromatic sp3 carbon.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or atom.GetIsAromatic() or atom.GetHybridization() != rdchem.HybridizationType.SP3:
                valid_ring = False
                break
        if not valid_ring:
            continue
        
        # Ensure that every bond between adjacent ring atoms is a single bond.
        bonds_ok = True
        for i in range(len(ring)):
            a1 = ring[i]
            a2 = ring[(i+1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(a1, a2)
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                bonds_ok = False
                break
        if not bonds_ok:
            continue
        
        # Check if any ring atom is substituted with a free carboxyl group.
        free_carboxyl_found = False
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 6 and is_free_carboxyl(nbr):
                    free_carboxyl_found = True
                    break
            if free_carboxyl_found:
                break
        
        if free_carboxyl_found:
            msg = "Found a cyclohexane ring with at least one attached free carboxyl substituent, consistent with a quinic acid core"
            return True, msg
    
    return False, "No cyclohexane ring with an attached free carboxyl substituent detected"
    
# Example test cases (uncomment to try):
# test_smiles = [
#     "O[C@H]1C[C@@](O)(C[C@H](O)[C@H]1O)C(O)=O",   # (+)-quinic acid (should be True)
#     "C1=CC(=C(C=C1)O)O"                             # Hydroquinone (should be False)
# ]
# for smi in test_smiles:
#     result, reason = is_quinic_acid(smi)
#     print(f"{smi}: {result} ({reason})")