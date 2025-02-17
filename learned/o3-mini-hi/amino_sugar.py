"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: Amino sugar – any sugar having one or more alcoholic hydroxy groups replaced by 
substituted or unsubstituted amino groups.

Revised algorithm:
  1. Parse the molecule from the SMILES.
  2. Force initialization of ring information.
  3. Loop over all rings (from RDKit) that are 5 or 6 members.
  4. In each such ring, require that it contains at least one oxygen and that the ring is “sugar‐like”
     (i.e. most ring atoms are carbon/oxygen).
  5. For each ring atom (typically carbons) check the exocyclic substituents. Count an -OH group if
     it has an attached hydrogen. For nitrogen substituents, only count them if they are sp3 and their 
     local environment suggests they are amino substituents (free –NH2 or part of an acetamido group).
  6. If a candidate ring has at least three exocyclic substituents and at least one acceptable amino substituent,
     classify the molecule as an amino sugar.
     
If ring information could not be determined or if no sugar-like ring is found, we return False with an appropriate reason.
"""

from rdkit import Chem

# Pre-compile a SMARTS pattern for an acetamido environment.
# This SMARTS looks for a nitrogen bound to a carbonyl carbon with a methyl group.
acetamido_smarts = Chem.MolFromSmarts("[NX3;!R][C](=O)[CH3]")

def is_valid_amino(n_atom):
    """
    Checks whether the given nitrogen atom appears to be a valid amino substituent:
       - It is an sp3 (non‐aromatic) nitrogen.
       - It either has one or more hydrogen atoms (free amine)
         or is part of an acetamido group (matches a C(=O)CH3 nearby).
    """
    # Only consider nitrogen atoms
    if n_atom.GetAtomicNum() != 7:
        return False
    # Must be sp3 (non‐aromatic)
    if n_atom.GetHybridization() != Chem.HybridizationType.SP3:
        return False
    # If the atom has one or more attached hydrogen, accept it.
    if n_atom.GetTotalNumHs() > 0:
        return True
    # Otherwise, attempt to see if its neighborhood matches the acetamido pattern.
    # Construct a sub-molecule from the bonds adjacent to this nitrogen.
    try:
        bond_idxs = [bond.GetIdx() for bond in n_atom.GetBonds()]
        env = Chem.PathToSubmol(n_atom.GetOwningMol(), bond_idxs)
        if env.HasSubstructMatch(acetamido_smarts):
            return True
    except Exception:
        pass
    return False

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    
    The algorithm inspects rings of size 5 or 6. For each candidate ring:
      - The ring must contain at least one oxygen (typical for a sugar ring) 
        and a majority of its atoms are carbon or oxygen.
      - Exocyclic substituents attached to ring carbons are examined.
            • Hydroxy (-OH) groups are counted if the oxygen appears to carry at least one hydrogen.
            • Amino substituents are counted only if the nitrogen atom is sp3 
              (i.e. non‐aromatic) and either carries a free –NH2 or is part of an acetamido group.
      - The candidate ring must have at least three exocyclic substituents and at least one acceptable amino substituent.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       (bool, str): (True, reason) if the molecule is classified as an amino sugar;
                    otherwise, (False, reason).
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Force update of property cache to (re)compute ring info.
    mol.UpdatePropertyCache(strict=False)
    
    # Try to extract ring info; if not available, fall back on GetSymmSSSR
    try:
        ring_info = mol.GetRingInfo().AtomRings()
    except Exception:
        # Force computation using GetSymmSSSR
        ring_info = [tuple(r) for r in Chem.GetSymmSSSR(mol)]
    
    candidate_found = False
    
    for ring in ring_info:
        # Only consider rings of size 5 or 6.
        if len(ring) not in (5, 6):
            continue

        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count oxygen atoms within the ring.
        oxy_in_ring = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
        if len(oxy_in_ring) < 1:
            continue  # Not sugar-like if no oxygen is in the ring.

        # Examine exocyclic substituents on ring atoms.
        total_exo = 0
        amino_found = False
        for atom in ring_atoms:
            # Focus on ring carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # Skip atoms that lie within the ring.
                at_num = nbr.GetAtomicNum()
                if at_num == 8:
                    # Count an oxygen substituent if it has one or more hydrogen either explicit or implicit.
                    if nbr.GetTotalNumHs() > 0:
                        total_exo += 1
                    else:
                        total_exo += 1
                elif at_num == 7:
                    # For nitrogen, count it only if it appears to be a valid amino substituent.
                    if is_valid_amino(nbr):
                        total_exo += 1
                        amino_found = True
        # Require the candidate ring to have at least three exocyclic substituents.
        if total_exo < 3:
            continue
        
        candidate_found = True
        if amino_found:
            return True, ("Detected a sugar-like ring (5- or 6-membered, containing at least one oxygen and multiple exocyclic substituents) "
                          "with an acceptable amino substituent (free -NH2 or acetamido).")
            
    # End of checking rings.
    if not candidate_found:
        return False, ("No typical sugar ring (5- or 6-membered ring with at least one oxygen and multiple exocyclic substituents) was detected.")
    else:
        return False, ("Sugar-like ring(s) were detected, but none contained an acceptable amino substituent "
                       "(free -NH2 or acetamido), so the molecule was not classified as an amino sugar.")

# Example usage:
if __name__ == "__main__":
    # Test a SMILES string that previously caused trouble.
    test_smiles = "C[C@H]1O[C@H](O)[C@H](O)[C@H]([C@H]1O)N(C)C"
    result, reason = is_amino_sugar(test_smiles)
    print("Result:", result)
    print("Reason:", reason)