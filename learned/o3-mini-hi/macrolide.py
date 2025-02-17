"""
Classifies: CHEBI:25106 macrolide
"""
#!/usr/bin/env python
"""
Classifies: macrolide
Definition: A macrolide is a macrocyclic lactone defined as a ring of 12 or more atoms 
            that contains an ester (lactone) function. In our algorithm the lactone is 
            identified by an oxygen–carbon bond inside a ring where the carbon is doubly 
            bonded to an oxygen that is not part of the ring. In addition, for a macrolide 
            (typically derived from polyketides) the ring should not contain nitrogen.
            
This revised algorithm:
  1. Parses the SMILES string.
  2. Retrieves all ring atom indices from the molecule.
  3. Loops over each ring of size >=12 and skips rings that contain any nitrogen atom.
  4. For each bond wholly contained in the ring, if the bond is between an oxygen and 
     a carbon then it inspects whether that carbon has a double bond to an oxygen that 
     is NOT in the ring (i.e. is exocyclic), which is taken as the lactone (ester) feature.
  5. If such a feature is found, it is classified as a macrolide.
  
If none of the rings show this feature, the molecule is not classified as a macrolide.
"""
from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines whether the molecule is a macrolide based on its SMILES string.
    A macrolide is defined as a macrocyclic lactone having a ring of 12 or more atoms 
    that features an ester bond (an oxygen–carbon bond where the carbon bears a double 
    bond to an oxygen not in the ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple (True/False, explanation) indicating classification result.
    """
    # Parse SMILES into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in the molecule"
    
    # Loop over each ring that is large enough.
    for ring in ring_info:
        if len(ring) < 12:
            continue  # Only consider rings of size 12 or more.
        ring_set = set(ring)
        # Exclude rings with nitrogen atoms (atomic num 7) because macrolide rings are typically 
        # polyketide-derived and predominantly contain carbon and oxygen.
        if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 7 for idx in ring_set):
            continue
        
        # Loop over bonds that are completely inside the ring.
        for bond in mol.GetBonds():
            idx1 = bond.GetBeginAtomIdx()
            idx2 = bond.GetEndAtomIdx()
            # Ensure both atoms are in the current ring.
            if idx1 not in ring_set or idx2 not in ring_set:
                continue
            
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # Identify candidate bond if one atom is oxygen and the other is carbon.
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8):
                # Let the carbon atom be the one with atomic number 6.
                carbon = atom2 if atom1.GetAtomicNum() == 8 else atom1
                # Now check if the carbon is carbonylated (has a double bond to an oxygen)
                # that is exocyclic (i.e. that oxygen is not in the ring).
                for nb_bond in carbon.GetBonds():
                    # Look for a double bond.
                    if nb_bond.GetBondTypeAsDouble() == 2:
                        other_atom = nb_bond.GetOtherAtom(carbon)
                        if other_atom.GetAtomicNum() == 8:
                            # Check that the carbonyl oxygen is not part of the ring.
                            if other_atom.GetIdx() not in ring_set:
                                reason = (f"Found macrocyclic lactone ring of size {len(ring)} with ester bond: "
                                          "oxygen in ring attached to carbonyl carbon (carbonyl oxygen exocyclic).")
                                return True, reason
    return False, "No macrocyclic lactone ring (ester embedded in a ring of 12 or more atoms) found"

# Example usage:
# test_smiles = "O1CCCCCCCCCC(OCCCCC1)=O"  # Adjust this example SMILES as needed.
# result, explanation = is_macrolide(test_smiles)
# print(result, explanation)