"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: macrolide
Definition: A macrolide is a macrocyclic lactone defined as a ring of 12 or more atoms 
            that features an ester (lactone) bond. Our algorithm identifies the lactone 
            feature by first finding a bond between an oxygen and a carbon where both atoms 
            are in a ring and where that carbon also bears a double bond to an oxygen that 
            is exocyclic (not part of the ring). In addition, we require that the macrocycle 
            (the ring containing the ester bond) does not contain any nitrogen atoms 
            (since macrolides are typically polyketide‐derived) and that the overall molecule 
            does not contain phosphorus (which would point to phospholipid‐like structures).
            
This version was revised after noticing that our previous approach both erroneously 
classified some non–macrolide structures (false positives) and missed some macrolides 
(false negatives). We therefore now use a targeted SMARTS query along with ring checks.
"""

from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines whether the molecule is a macrolide based on its SMILES string.
    A macrolide is defined as a macrocyclic lactone having a ring of 12 or more atoms 
    that features an ester bond: an oxygen in the ring connected to a carbon that bears 
    a double bond to an exocyclic oxygen. Also, the ring should not contain nitrogen and 
    the molecule should not contain phosphorus.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple (True/False, explanation) indicating the classification result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that contain phosphorus (atomic number 15),
    # which are not typical of polyketide-derived macrolides.
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Molecule contains phosphorus; likely not a macrolide (polyketide-derived) structure."
    
    # Define a SMARTS pattern for the lactone bond:
    #   [O;r] - an oxygen that is a member of a ring, bonded to
    #   [C;r](=[O;!r]) - a carbon (in a ring) that is double-bonded to an oxygen that is NOT in a ring.
    lactone_smarts = "[O;r]-[C;r](=[O;!r])"
    lactone_pat = Chem.MolFromSmarts(lactone_smarts)
    if lactone_pat is None:
        return False, "Error building SMARTS pattern"

    # Find all matches of lactone substructure.
    matches = mol.GetSubstructMatches(lactone_pat)
    if not matches:
        return False, "No lactone substructure (ester bond with exocyclic carbonyl oxygen) found."
    
    # Retrieve ring information from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Loop over each lactone substructure match.
    # The match gives a tuple of atom indices: (ring oxygen, ring carbon, carbonyl oxygen).
    for match in matches:
        ring_o, ring_c, exo_o = match
        # Now look for a ring (from RDKit's ring perception) that contains both ring_o and ring_c.
        for ring in ring_info:
            if len(ring) < 12:
                continue  # Only interested in macrocycles (size 12 or more).
            if (ring_o in ring) and (ring_c in ring):
                # For our purposes the ring itself must not contain any nitrogen.
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring):
                    continue
                reason = (f"Found macrocyclic lactone ring of size {len(ring)} with ester bond: "
                          "oxygen in ring attached to carbonyl carbon (carbonyl oxygen exocyclic).")
                return True, reason

    return False, "No macrocyclic lactone ring (ester embedded in a ring of 12 or more atoms) found"

# Example usage:
# test_smiles = "O1CCCCCCCCCC(OCCCCC1)=O"  # A tentative macrolide example.
# result, explanation = is_macrolide(test_smiles)
# print(result, explanation)