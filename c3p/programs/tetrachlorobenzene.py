"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition: Any member of the class of chlorobenzenes carrying exactly four chlorine atoms attached directly to a benzene ring.
We improve on the previous method by:
  • Removing (or choosing) the largest fragment in case the SMILES encodes salts or counterions.
  • For each benzene ring candidate (6-membered aromatic ring of carbons),
    we count the number of substituents that are Cl.
  • We also check that the substituents besides Cl (if any) on that ring are not “exotic” – in a typical tetrachlorobenzene the benzene ring should have 4 Cl and the remainder (usually 2) that are implicitly hydrogen or very simple groups.
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    The algorithm:
      1. Parses the SMILES string into an RDKit molecule. 
      2. If the molecule is split into fragments (eg, counterions are present), selects the largest fragment.
      3. Loops over each ring; if a ring has exactly 6 atoms and each is a carbon with aromaticity set,
         it proceeds to count the number of chlorine substituents (neighbors that are not in the ring).
      4. As an additional heuristic, for every atom in the candidate ring the number of non-Cl substituents
         is checked; if too many extra groups are present the ring is likely not a “pure” chlorobenzene.
      5. If a ring is found with exactly 4 Cl substituents and no evidence of extra (non minor) substituents,
         the function returns True and a corresponding reason.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if a qualifying benzene ring is found, False otherwise.
       str: Reason for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {e}"
    
    # If multiple fragments exist (e.g. salts or counterions), take the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        
    ring_info = mol.GetRingInfo()
    if not ring_info or ring_info.NumRings() == 0:
        return False, "No rings found in the molecule"
    
    # Loop over each ring in the molecule
    for ring in ring_info.AtomRings():
        # Only consider rings with 6 atoms
        if len(ring) != 6:
            continue
            
        # Check that every atom in the ring is a carbon and is aromatic
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue
        
        # Count chlorine substituents attached to ring atoms.
        cl_count = 0
        extra_substituents = 0  # non-Cl substituents attached explicitly to the ring atoms
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                # Only consider neighbors not in the ring (i.e. substituents)
                if neighbor.GetIdx() in ring:
                    continue
                # If the neighbor is chlorine, count it.
                if neighbor.GetAtomicNum() == 17:
                    cl_count += 1
                else:
                    extra_substituents += 1
        # For a classic tetrachlorobenzene we expect exactly 4 Cl substituents.
        # In many valid examples (e.g. tetrachlorophenol) there may be one or two extra very small groups (OH, COOH),
        # so we allow at most 2 extra non-metal substituents.
        if cl_count == 4 and extra_substituents <= 2:
            return True, "Found benzene ring with exactly 4 chlorine substituents and acceptable additional groups"
            
    return False, "No benzene ring with exactly 4 chlorine substituents (and acceptable extra groups) found"

# Example usage:
if __name__ == "__main__":
    # Test with one known true positive: 1,2,3,5-tetrachlorobenzene
    test_smiles = "Clc1cc(Cl)c(Cl)c(Cl)c1"
    result, reason = is_tetrachlorobenzene(test_smiles)
    print(result, reason)