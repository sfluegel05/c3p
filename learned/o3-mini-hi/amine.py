"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
Definition: A compound formally derived from ammonia by replacing one, two or three hydrogen atoms 
by hydrocarbyl groups. Primary (–NH2), secondary (–NH–) or tertiary (–N–) amines are valid. 
A nitrogen directly bonded to a carbonyl (e.g. in an amide) or that is part of a heteroaromatic ring 
(e.g. pyridine) is excluded.
"""

from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule contains a qualifying amine functional group
    based on its SMILES string. The algorithm iterates over nitrogen atoms and checks:
      - The nitrogen atom is not quaternary.
      - It is not part of a heteroaromatic ring (e.g. pyridine) but it is acceptable if it is exocyclic.
      - It is not directly attached to a carbonyl carbon (i.e. a typical amide bond).
      - It has a hydrogen count consistent with a primary (2 H), secondary (1 H) 
        or tertiary (0 H) derivative of ammonia.
        
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if at least one qualifying free-amine group is found, False otherwise.
        str: An explanation (which type of amine was detected or why not).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate over all atoms looking for nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # not a nitrogen
        # Skip if formal charge is positive (e.g. quaternary)
        if atom.GetFormalCharge() > 0:
            continue
        # In many heterocyclic systems (eg pyridine) N is aromatic and in-ring.
        # We assume that a nitrogen atom that is in an aromatic ring is not a free amine.
        if atom.GetIsAromatic() and atom.IsInRing():
            continue

        # Get the total number of hydrogens attached to this nitrogen.
        nH = atom.GetTotalNumHs()
        # We also want to ignore cases where the nitrogen appears “over‐coordinated”
        if atom.GetDegree() > 3:
            continue

        # Check if the nitrogen is directly attached to a carbonyl carbon.
        attached_to_carbonyl = False
        for neighbor in atom.GetNeighbors():
            # Only consider carbon neighbors
            if neighbor.GetAtomicNum() == 6:
                for bond in neighbor.GetBonds():
                    # Look for a double bond from the carbon neighbor to an oxygen.
                    if bond.IsDouble():
                        other = bond.GetOtherAtom(neighbor)
                        if other.GetAtomicNum() == 8:
                            attached_to_carbonyl = True
                            break
                if attached_to_carbonyl:
                    break

        # If the N is directly bound to a carbonyl, then we do not count it.
        if attached_to_carbonyl:
            continue

        # At this point, we have a nitrogen not in an amide or heteroaromatic ring.
        # Now classify the type based on the number of hydrogen atoms.
        if nH == 2:
            return True, "Found a primary amine (–NH2) group not bound to a carbonyl."
        elif nH == 1:
            return True, "Found a secondary amine (–NH–) group not bound to a carbonyl."
        elif nH == 0:
            # Tertiary amine: note that sometimes the nitrogen in aromatic amines might have 0 H.
            # We already skipped in-ring aromatic N so here we assume it is a tertiary alkyl amine.
            return True, "Found a tertiary amine (–N–) group not bound to a carbonyl."
        # If unexpected hydrogen count, continue checking other nitrogens.
    
    # If we see no qualifying nitrogen, return false.
    return False, "No qualifying free amine functional group found."

# Example usage:
# result, reason = is_amine("CNc1ccccc1")
# print(result, reason)