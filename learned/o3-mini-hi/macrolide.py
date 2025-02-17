"""
Classifies: CHEBI:25106 macrolide
"""
"""
Classifies: macrolide
Definition: A macrolide is a macrocyclic lactone defined as a ring of 12 or more atoms that contains an ester bond 
             (i.e. an oxygen–connected carbonyl group) as part of the ring. Such molecules are typically derived from 
             polyketide biosynthesis.
             
This improved algorithm:
  1. Parses the SMILES string.
  2. Retrieves all rings (using RDKit’s ring information).
  3. For each ring of size 12 or more, inspects bonds among atoms in the ring.
     It looks for a bond connecting an oxygen atom (the “ester oxygen” in the ring) with a carbon atom, for which 
     the carbon atom has at least one double bond (i.e. bond order 2) to another oxygen (the carbonyl oxygen).
  4. If such a bond is found in any large ring, the molecule is classified as a macrolide.
  
Note: This does not attempt to determine polyketide origin and uses only the structure of a lactone ring.
"""
from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines whether a molecule is a macrolide based on its SMILES string.
    
    A macrolide is defined as a macrocyclic lactone having a ring of 12 or more atoms that features an ester bond.
    The ester (lactone) bond is identified here by looking for a bond between an oxygen (in the ring) and a carbon 
    that is double-bonded to an oxygen (the carbonyl). 
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple containing a boolean classification and a textual explanation.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get the ring information for the molecule.
    ring_atom_tuples = mol.GetRingInfo().AtomRings()  # Each tuple contains atom indices for a ring.
    if not ring_atom_tuples:
        return False, "No rings found in the molecule"
    
    # Loop over rings that are large enough (>=12 atoms).
    for ring in ring_atom_tuples:
        if len(ring) < 12:
            continue  # Skip rings that are too small.
            
        # Create a set for faster membership checking.
        ring_set = set(ring)
        
        # For each bond in the molecule that is entirely inside the ring, check for an ester bond feature.
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            # Check that both atoms of the bond are in the current ring.
            if idx1 not in ring_set or idx2 not in ring_set:
                continue
            
            # Check if one of the two is an oxygen and the other is a carbon.
            # (This bond is a candidate for being the ring-embedded "ester" linkage.)
            if (a1.GetAtomicNum() == 8 and a2.GetAtomicNum() == 6) or (a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 8):
                # Identify the carbon atom.
                carbon_atom = a2 if a1.GetAtomicNum() == 8 else a1
                
                # Now check if this carbon is carbonylated (has a double bond to an oxygen)
                # Loop over its bonds.
                for nb_bond in carbon_atom.GetBonds():
                    # Ensure the bond order is double.
                    if nb_bond.GetBondTypeAsDouble() == 2:
                        other_atom = nb_bond.GetOtherAtom(carbon_atom)
                        if other_atom.GetAtomicNum() == 8:
                            # Found a carbonyl double bond.
                            reason = f"Found macrocyclic lactone ring of size {len(ring)} with ester bond: oxygen (in ring) attached to carbonyl carbon."
                            return True, reason

    return False, "No macrocyclic lactone ring (ester embedded in a ring of 12 or more atoms) found"

# Example usage (for testing):
# test_smiles = "O1CCCCCCCCCC(OCCCCC1)=O"  # Adjust example as needed
# result, explanation = is_macrolide(test_smiles)
# print(result, explanation)