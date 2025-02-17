"""
Classifies: CHEBI:36092 clavulone
"""
"""
Classifies: Clavulone – a class of esterified prostanoids obtained from marine corals.
Heuristic criteria (improved):
  1. The molecule must contain at least one 5-membered ring that has:
       • at least one exocyclic carbonyl (oxygen double-bonded to a ring atom that is not part of the ring)
         and among these, at least one is conjugated with an internal double bond,
       • exactly one double bond between ring atoms,
       • and the ring is not fused with any other ring (i.e. does not share >1 atom with any other ring).
  2. The molecule must contain at least one ester substituent
       (defined by the substructure “[#6]-[OX2]-[CX3](=O)[#6]”) that is within a graph distance ≤2 of
       at least one atom in the candidate ring.
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_clavulone(smiles: str):
    """
    Determines if a molecule qualifies as a clavulone (esterified prostanoid) based on its SMILES string.
    
    The criteria (heuristic) are:
      1. The molecule contains at least one 5-membered ring that:
           - has at least one exocyclic carbonyl (i.e. an oxygen atom double-bonded to a ring atom while
             not being part of the ring),
           - has exactly one (internal) double bond between ring atoms,
           - and at least one of the exocyclic carbonyl-bearing ring atoms is part of that internal double bond,
           - with the ring not being fused with another ring.
      2. The molecule contains at least one ester substituent (SMARTS: “[#6]-[OX2]-[CX3](=O)[#6]”) 
         that is in close proximity (graph distance ≤2 bonds) to any atom in the candidate ring.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule meets the clavulone heuristic, False otherwise.
        str : Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    candidate_ring = None  # will hold the indices for the candidate 5-membered ring
    
    # Loop over 5-membered rings
    for ring in all_rings:
        if len(ring) != 5:
            continue
        
        ring_atom_set = set(ring)
        
        # 1a. Identify exocyclic carbonyl groups:
        # For each atom in the ring, check if it has a neighbor (not in the ring) that is oxygen
        # attached by a double bond.
        carbonyl_atoms = []  # store ring atom indices that bear an exocyclic carbonyl
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring_atom_set:
                    continue
                # Check for oxygen attached with a double bond
                if nbr.GetAtomicNum() == 8:
                    bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
                    if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                        carbonyl_atoms.append(atom_idx)
        if len(carbonyl_atoms) < 1:
            # no exocyclic carbonyl found in this ring
            continue
            
        # 1b. Count internal (ring–ring) double bonds.
        internal_double_bonds = []
        # Loop over all bonds in the mol; only consider bonds connecting atoms in the ring.
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_atom_set and a2 in ring_atom_set:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    internal_double_bonds.append((a1, a2))
        if len(internal_double_bonds) != 1:
            continue
        
        # 1c. Check that one of the ring atoms bearing an exocyclic carbonyl is part of the internal double bond.
        db_atoms = set(internal_double_bonds[0])
        conjugated = any(atom in db_atoms for atom in carbonyl_atoms)
        if not conjugated:
            continue
        
        # 1d. Check that the candidate ring is not fused with any other ring.
        fused = False
        for other_ring in all_rings:
            if other_ring == ring:
                continue
            # if the intersection of atoms is more than 1, the ring is fused
            if len(set(ring).intersection(set(other_ring))) > 1:
                fused = True
                break
        if fused:
            continue
        
        # If all criteria are met, we have found our candidate ring.
        candidate_ring = list(ring)
        break
    
    if candidate_ring is None:
        msg = ("No suitable cyclopentenone-type ring found: need a 5-membered ring with at least one exocyclic "
               "carbonyl (conjugated to the single internal double bond) and not fused with another ring.")
        return False, msg
    
    # 2. Look for at least one ester substituent.
    # Define the ester SMARTS pattern.
    ester_smarts = "[#6]-[OX2]-[CX3](=O)[#6]"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    if ester_query is None:
        return False, "Error generating ester SMARTS pattern"
        
    ester_matches = mol.GetSubstructMatches(ester_query)
    ester_found = False
    # Check if any atom in an ester match is within graph distance 2 (i.e. 2 bonds) of any candidate ring atom.
    for match in ester_matches:
        for match_atom in match:  # match is a tuple of atom indices
            for ring_atom in candidate_ring:
                path = rdmolops.GetShortestPath(mol, ring_atom, match_atom)
                if len(path) - 1 <= 2:
                    ester_found = True
                    break
            if ester_found:
                break
        if ester_found:
            break
    
    if not ester_found:
        msg = ("Candidate cyclopentenone-type ring found but no ester substituent found in close proximity "
               "(graph distance ≤2) to the ring.")
        return False, msg
    
    msg = ("Contains a cyclopentenone-type 5-membered ring with an exocyclic carbonyl conjugated to its internal double bond "
           "plus an ester substituent in close proximity – consistent with clavulone.")
    return True, msg

# Example usage (for testing):
if __name__ == "__main__":
    # Testing with one of the clavulone examples: punaglandin 2
    smiles_example = ("ClC=1C(=O)[C@@]([C@@](O)(C/C=C\\CCCCC)C1)"
                      "([C@@H](OC(=O)C)[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O)[H]")
    result, reason = is_clavulone(smiles_example)
    print(result, reason)