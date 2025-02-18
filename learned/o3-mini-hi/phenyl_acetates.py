"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl Acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid 
with the hydroxy group of any phenol.
We require that the ester fragment –O–C(=O)CH3 is directly attached to an aromatic carbon,
and importantly, that the ester bond is acyclic (i.e. not part of a lactone).
"""

from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines whether a molecule is a phenyl acetate based on its SMILES string.
    
    A phenyl acetate contains an aromatic (phenolic) oxygen that has been acetylated,
    forming an ester bond: R-c-O-C(=O)CH3.
    In this algorithm we:
      1. Look for the substructure defined by the SMARTS "[c:1][O:2][C:3](=O)[CH3]".
      2. For each match, check that the bond between the ester oxygen and the carbonyl carbon
         is not part of a ring (to avoid cyclic esters like lactones).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as a phenyl acetate, False otherwise.
        str: Explanation for the classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for an acetyl group directly attached to an aromatic (phenolic) oxygen.
    # [c:1] represents an aromatic carbon atom (part of a phenol), and [O:2] is the oxygen,
    # [C:3](=O) is the carbonyl, and [CH3] is the acetyl methyl group.
    acetate_pattern = Chem.MolFromSmarts("[c:1][O:2][C:3](=O)[CH3]")
    if acetate_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    matches = mol.GetSubstructMatches(acetate_pattern)
    if not matches:
        return False, "No phenyl acetate substructure ([c][O][C](=O)C) found in the molecule."
    
    # Get the ring information from the molecule; this will help us
    # ensure that the ester bond is not part of a ring.
    ring_info = mol.GetRingInfo().AtomRings()
    
    # For each match, verify that the ester oxygen (atom index  match[1])
    # and the carbonyl carbon (atom index match[2]) are not in a same ring.
    valid_match_found = False
    for match in matches:
        # match is a tuple (idx_c, idx_O, idx_Ccarbonyl, idx_CH3)
        idx_oxygen = match[1]
        idx_carbonyl = match[2]
        
        in_same_ring = False
        for ring in ring_info:
            if idx_oxygen in ring and idx_carbonyl in ring:
                in_same_ring = True
                break
        
        if in_same_ring:
            # This match likely corresponds to a cyclic ester (lactone) rather than a phenyl acetate.
            continue
        
        # Optionally, we might check that the aromatic carbon (match[0]) is indeed aromatic
        # and has at least one hydrogen; here RDKit already assigns aromaticity.
        atom_aromatic = mol.GetAtomWithIdx(match[0])
        if not atom_aromatic.GetIsAromatic():
            continue
            
        valid_match_found = True
        break

    if valid_match_found:
        return True, "Molecule contains an acyclic phenyl acetate moiety (acetylated phenolic oxygen)."
    else:
        return False, "Found acetyl ester fragment but it is embedded in a ring system (likely a lactone), not a phenyl acetate."

# Example usage:
if __name__ == "__main__":
    # Test with phenyl acetate itself
    test_smiles = "CC(=O)Oc1ccccc1"
    result, reason = is_phenyl_acetates(test_smiles)
    print(result, reason)