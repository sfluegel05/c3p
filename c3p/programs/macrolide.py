"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: macrolide
Definition: A macrolide is defined as a macrocyclic lactone having a ring of twelve or more
            atoms featuring an ester bond. More specifically, an ester oxygen (within the ring)
            is connected to a carbon that is double-bonded to an exocyclic oxygen. The macrocycle
            should not contain nitrogen atoms and the overall molecule should not contain phosphorus.
            
Improvements:
  • Only the ester oxygen is required to be in a ring; then we search through the molecule’s
    ring information to find a macrocycle (≥12 atoms) that contains the ester oxygen and the
    attached carbon while excluding the carbonyl oxygen.
  • The identified ring is further filtered to contain no nitrogen atoms and be predominantly
    composed of carbon atoms (at least 60%).
  • Molecules containing any phosphorus atoms are rejected upfront.
"""

from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines whether the molecule is a macrolide based on its SMILES string.
    A macrolide is defined as a macrocyclic lactone (ring of >=12 atoms) that has an ester bond 
    with an oxygen atom (ester oxygen) that is in a ring and bonded to a carbon that is double-bonded 
    to an exocyclic oxygen (carbonyl oxygen). Additionally, the ring should not contain nitrogen atoms
    and the molecule must not contain any phosphorus.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): A tuple with a boolean (True if macrolide, else False) and a reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules containing phosphorus (atomic number 15)
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely not a macrolide (polyketide-derived) structure."
    
    # Build a SMARTS pattern for a lactone: an ester oxygen in a ring bonded to a carbon that has a double bond to an oxygen.
    # Note: We do not require that the carbonyl oxygen be in the same ring.
    lactone_smarts = "[O;r]-[C](=[O])"
    lactone_pat = Chem.MolFromSmarts(lactone_smarts)
    if lactone_pat is None:
        return False, "Error building SMARTS pattern"
    
    # Look for any lactone substructure match.
    matches = mol.GetSubstructMatches(lactone_pat)
    if not matches:
        return False, "No lactone substructure (ester bond with exocyclic carbonyl) found."
    
    # Get all ring definitions for the molecule.
    rings = mol.GetRingInfo().AtomRings()
    if not rings:
        return False, "No rings found in the molecule."
    
    # Iterate over all lactone matches.
    # Each match is a tuple of indices: (ester oxygen index, adjacent carbon index, carbonyl oxygen index)
    for match in matches:
        est_o_idx, c_idx, cO_idx = match
        
        # For an acceptable macrolide, we need to find at least one ring (with 12 or more atoms)
        # that contains both the ester oxygen and its bonded carbon, but does NOT contain the carbonyl oxygen.
        for ring in rings:
            if len(ring) < 12:
                continue  # Not a macrocyclic ring.
            if (est_o_idx in ring) and (c_idx in ring) and (cO_idx not in ring):
                # Reject if the ring contains any nitrogen (atomic number 7)
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
                    continue
                # Ensure that at least 60% of the atoms in the ring are carbons (atomic number 6)
                n_carbons = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
                if (n_carbons / len(ring)) < 0.60:
                    continue  # Too many heteroatoms; not typical for a polyketide macrocycle.
                
                reason = (f"Found macrocyclic lactone ring of size {len(ring)} with ester bond: "
                          f"ester oxygen (atom {est_o_idx}) and carbon (atom {c_idx}) are in the ring, "
                          f"with exocyclic carbonyl oxygen (atom {cO_idx}); ring composition is {n_carbons}/{len(ring)} carbons.")
                return True, reason

    return False, "No macrocyclic lactone ring (ester embedded in a ring of 12 or more atoms) found."

# Example usage:
# test_smiles = "O1CCCCCCCCCC(OCCCCC1)=O"  # tentative macrolide candidate.
# result, explanation = is_macrolide(test_smiles)
# print(result, explanation)