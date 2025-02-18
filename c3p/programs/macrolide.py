"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: macrolide
Definition: A macrolide is defined as a macrocyclic lactone having a ring of twelve or more 
            atoms that features an ester bond – that is, an oxygen in the ring connected to a 
            carbon (which in turn is double-bonded to an exocyclic oxygen). In addition, the 
            macrocycle should not contain nitrogen atoms (since macrolides are typically 
            polyketide‐derived) and the overall molecule should not contain phosphorus.
            
Reasoning improvements:
  • Instead of enforcing that both oxygen and carbon be “in ring” via SMARTS, we require only 
    that the ester oxygen is in a ring. Then for every ester match we search the ring info for a 
    ring (size ≥12 atoms) that contains the ester oxygen and its attached carbon.
  • In addition, we check that the double‐bonded (carbonyl) oxygen is not part of that same ring.
  • To filter out false positives from e.g. glycosides, we further require that the chosen macrocyclic 
    ring has at least 60% carbon atoms.
  • Molecules containing phosphorus are rejected.
"""

from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines whether the molecule is a macrolide based on its SMILES string.
    A macrolide is defined as a macrocyclic lactone (ring of >=12 atoms) that has an ester bond 
    with an oxygen atom that is in the ring and a carbon atom that is double-bonded to an oxygen 
    which is exocyclic. Also, the ring should not contain nitrogen and the molecule should not 
    contain phosphorus.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): A tuple with a boolean (True if macrolide, else False) and an explanation.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that contain phosphorus (atomic number 15).
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely not a macrolide (polyketide-derived) structure."
    
    # We start by looking for lactone substructures.
    # Our relaxed SMARTS requires:
    #  - An oxygen (the ester oxygen) that is a member of a ring [O;r]
    #  - Bonded to a carbon (not forcing ring membership on C) which is double-bonded to an oxygen.
    #    We do not restrict the second oxygen to be in a ring (expected exocyclic carbonyl).
    lactone_smarts = "[O;r]-[C](=[O])"
    lactone_pat = Chem.MolFromSmarts(lactone_smarts)
    if lactone_pat is None:
        return False, "Error building SMARTS pattern"
    
    matches = mol.GetSubstructMatches(lactone_pat)
    if not matches:
        return False, "No lactone substructure (ester bond with exocyclic carbonyl) found."
    
    # Get ring information from RDKit.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in the molecule."

    # We will iterate over all the lactone matches.
    # Each match returns a tuple: (ester oxygen index, carbon index, carbonyl oxygen index)
    for match in matches:
        est_o, c_atom, cO = match

        # Now attempt to find at least one ring (macrocycle) that:
        #   1. Has size >= 12.
        #   2. Contains both the ester oxygen and the adjacent carbon.
        #   3. Does NOT contain the carbonyl oxygen atom.
        #   4. Does not have any nitrogen atoms.
        #   5. Has at least 60% carbons.
        for ring in rings:
            if len(ring) < 12:
                continue  # Only interested in macrocycles.
            if (est_o in ring) and (c_atom in ring) and (cO not in ring):
                # Reject if any atom in the ring is nitrogen.
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
                    continue
                # Count carbon atoms in the ring.
                n_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                if (n_carbons / len(ring)) < 0.60:
                    continue  # Too many heteroatoms – not a typical polyketide macrocycle.
                reason = (f"Found macrocyclic lactone ring of size {len(ring)} with ester bond: "
                          f"ester oxygen (atom {est_o}) and carbon (atom {c_atom}) are in the ring "
                          f"with exocyclic carbonyl oxygen (atom {cO}); ring is {n_carbons}/{len(ring)} carbons.")
                return True, reason

    return False, "No macrocyclic lactone ring (ester embedded in a ring of 12 or more atoms) found."

# Example usage:
# test_smiles = "O1CCCCCCCCCC(OCCCCC1)=O"  # tentative macrolide example.
# result, explanation = is_macrolide(test_smiles)
# print(result, explanation)