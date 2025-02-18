"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition: Any member of the class of chlorobenzenes carrying exactly four chlorine atoms attached directly 
to a benzene ring. Valid molecules may have an extra substituent if it is either a simple group (e.g. –OH, –CN) or 
a benzene ring linkage as in biphenyls. Any extra substituent that is “exotic” (for example, one carrying a nonzero 
formal charge) disfavors the classification.
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES. If multiple fragments are present (salts, counterions), select the largest.
      2. Loop through each ring that is exactly 6 atoms:
         a. Require that every atom in the ring is an aromatic carbon (benzene candidate).
         b. For every ring atom, examine each neighbor (i.e. substituent not in the ring).
            - If the neighbor is chlorine (atomic number 17), count it.
            - Otherwise, assess whether the substituent is “allowed”. In our heuristic, an allowed substituent is:
                • An atom with zero formal charge and that is very simple. This “simplicity” is evaluated as:
                     - If it is oxygen and its degree is 1 (eg –OH) 
                     - If it is a carbon that is part of an aromatic ring (as in biphenyl linkages)
                     - Or if it is a carbon whose degree is 1 (as in a cyano group, where that C is triple-bonded to N)
                • Any substituent failing these tests is counted as an extra group.
         c. We finally require that the candidate benzene ring has exactly 4 chlorine substituents
            and no “non‐allowed” extra substituents.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if a qualifying benzene ring is found, False otherwise.
       str: Reason for the decision.
    """
    # Parse SMILES and try sanitizing.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {e}"
    
    # If multiple fragments exist (counterions etc.), select the largest fragment
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        
    ring_info = mol.GetRingInfo()
    if not ring_info or ring_info.NumRings() == 0:
        return False, "No rings found in the molecule"
    
    # Loop over each ring in the molecule
    for ring in ring_info.AtomRings():
        # Only consider rings with 6 atoms.
        if len(ring) != 6:
            continue
        
        # Check that every atom in the ring is a carbon and aromatic.
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue
        
        # For this benzene candidate, count Cl and examine non-Cl substituents.
        cl_count = 0
        non_allowed_extras = 0
        
        # Loop over each atom of the ring.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Consider each neighbor not in the ring.
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                # If neighbor is chlorine, count it.
                if neighbor.GetAtomicNum() == 17:
                    cl_count += 1
                else:
                    # If the substituent carries a nonzero formal charge, disallow.
                    if neighbor.GetFormalCharge() != 0:
                        non_allowed_extras += 1
                        continue
                    # Allow a simple oxygen atom (e.g. –OH)
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        continue
                    # Allow aromatic carbons (e.g. in biphenyl links)
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                        continue
                    # Allow a carbon substituent if it is very simple (degree==1, e.g. in a cyano group)
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1:
                        continue
                    # Otherwise, count as extra.
                    non_allowed_extras += 1
        
        # Check whether this ring qualifies as a tetrachlorobenzene.
        # It must have exactly 4 Cl substituents and no non-allowed extra substituents.
        if cl_count == 4 and non_allowed_extras == 0:
            return True, "Found benzene ring with exactly 4 chlorine substituents and acceptable additional groups"
            
    return False, "No benzene ring with exactly 4 chlorine substituents (and acceptable extra groups) found"

# Example usage:
if __name__ == "__main__":
    # One known true positive: 1,2,3,5-tetrachlorobenzene SMILES.
    test_smiles = "Clc1cc(Cl)c(Cl)c(Cl)c1"
    result, reason = is_tetrachlorobenzene(test_smiles)
    print(result, reason)