"""
Classifies: CHEBI:15734 primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Primary alcohol
Definition: A primary alcohol is a compound in which a hydroxyl group (-OH) is attached 
           to a saturated (sp3) carbon that is primary – having either no or only one other carbon attached.
           
In other words, the –OH must be bonded to a carbon that is either:
 • methanol-type: CH3–OH  (the carbon is attached only to hydrogens)
 • CH2OH-type: R–CH2–OH  (the carbon has exactly one other carbon attached)
           
This improved algorithm:
  1. Parses a SMILES string and adds explicit hydrogens.
  2. Loops over oxygen atoms that bear at least one hydrogen.
  3. For each –OH candidate, ensures that:
       - The oxygen is connected to exactly one heavy atom (the candidate carbon).
       - The candidate carbon is sp³ hybridized.
       - When ignoring the –OH oxygen, the candidate carbon is attached to at most one other carbon.
  4. Returns True with an appropriate reason if a primary alcohol group is found.
  5. Returns False otherwise.
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule contains a primary alcohol group based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple; first element is True if a primary alcohol is detected,
                     otherwise False; second element explains the result.
    """
    # Parse the SMILES and add explicit hydrogens, so the -OH hydrogens are visible
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Loop over all atoms to find candidate -OH groups.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        
        # Check that this oxygen has at least one hydrogen attached.
        if atom.GetTotalNumHs() < 1:
            continue
        
        # The -OH oxygen should have only one heavy (non-hydrogen) neighbor.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue
        
        # The heavy neighbor candidate should be a carbon.
        cand = heavy_neighbors[0]
        if cand.GetSymbol() != "C":
            continue
        
        # Ensure that the carbon is sp3 hybridized (i.e. saturated).
        if cand.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Count how many carbon neighbors (heavy atoms) the candidate carbon has,
        # ignoring the current -OH oxygen.
        c_neighbors = [nbr for nbr in cand.GetNeighbors() if nbr.GetSymbol() == "C" and nbr.GetIdx() != atom.GetIdx()]
        
        # For a primary alcohol, the candidate carbon should have at most one other carbon neighbor.
        if len(c_neighbors) > 1:
            continue
        
        # We now have a primary carbon bearing an OH.
        # Optionally, we can inspect the total hydrogen count on the candidate.
        # For methanol (CH3OH) we expect 3 hydrogens and for CH2OH we expect 2.
        h_count = cand.GetTotalNumHs()
        if len(c_neighbors) == 0:
            return True, "Found primary alcohol group (methanol-type: CH3OH)"
        elif len(c_neighbors) == 1 and h_count >= 1:
            return True, "Found primary alcohol group (CH2OH-type)"
        else:
            # Even if the hydrogen count is unexpected, we still accept the candidate.
            return True, "Found primary alcohol group (primary C–OH detected)"
            
    return False, "No primary alcohol group found"


# Optional: Testing the function on a subset of example SMILES
if __name__ == "__main__":
    test_examples = [
        ("(2S)-2-aminobutan-1-ol(1+)", "CC[C@@H](CO)[NH3+]"),
        ("Semiplenamide A", "O=C(NCCO)/C(=C/CC/C=C/CCCCCCCCCCCCC)/C"),
        ("3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,10-heptadecafluoro-1-decanol", "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
        ("furfuryl alcohol", "OCc1ccco1"),
        ("Ethanolamine", "NCCO"),
        ("docosan-1-ol", "OCCCCCCCCCCCCCCCCCCCCCO")
    ]
    
    for name, smi in test_examples:
        result, reason = is_primary_alcohol(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")