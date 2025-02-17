"""
Classifies: CHEBI:15734 primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Primary Alcohol

Definition: A primary alcohol is a compound in which a hydroxyl group (-OH) is attached 
to a saturated (sp3) carbon that is primary – i.e. either:
  • methanol-type: the carbon has three H’s (CH3OH), or 
  • CH2OH-type: the carbon is attached to one (and only one) other carbon and bears two H’s.
  
The algorithm:
  1. Parses the SMILES string and adds explicit hydrogens.
  2. Iterates over all oxygen atoms and selects those that have at least one explicit hydrogen,
     and only one heavy (non-H) neighbor.
  3. For the candidate oxygen’s heavy neighbor (which must be carbon):
       a. Check that it is sp3-hybridized.
       b. Count its explicit hydrogen atoms.
       c. Count how many other carbons it is attached to (ignoring the –OH oxygen).
  4. If the candidate carbon has no carbon neighbors and exactly 3 hydrogens, it is a CH3OH group.
     Or if it has exactly one carbon neighbor and exactly 2 hydrogens, it is a CH2OH group.
  5. Return True if any such group is found; otherwise, return False.
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
    # Parse SMILES and add explicit hydrogens so that OH hydrogens are visible.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Loop over atoms to find an oxygen candidate for a primary –OH group.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        
        # Count explicit hydrogens attached to this oxygen.
        oh_hydrogens = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "H"]
        if len(oh_hydrogens) < 1:
            continue  # This oxygen does not bear any –OH hydrogen.
        
        # The oxygen in an alcohol should be connected to exactly one heavy atom.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue
        
        # The heavy neighbor should be a carbon.
        cand_carbon = heavy_neighbors[0]
        if cand_carbon.GetSymbol() != "C":
            continue
        
        # Check that the candidate carbon is sp3 hybridized (i.e. saturated).
        if cand_carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Count the number of explicit hydrogens on the candidate carbon.
        # (We use explicit atoms; note that Chem.AddHs already made them explicit.)
        c_hcount = len([nbr for nbr in cand_carbon.GetNeighbors() if nbr.GetSymbol() == "H"])
        
        # Count the number of other carbon neighbors (ignoring the current –OH oxygen).
        other_carbons = [nbr for nbr in cand_carbon.GetNeighbors() 
                         if nbr.GetSymbol() == "C" and nbr.GetIdx() != atom.GetIdx()]
        
        # For a primary alcohol the candidate carbon must be
        # either methyl (CH3OH) or methylene (CH2OH) attached to one carbon.
        if len(other_carbons) == 0 and c_hcount == 3:
            return True, "Found primary alcohol group (methanol-type: CH3OH)"
        elif len(other_carbons) == 1 and c_hcount == 2:
            return True, "Found primary alcohol group (CH2OH-type)"
        # Otherwise, not a candidate for a primary alcohol.
        # (Even if the candidate carbon bears unexpected hydrogen counts, we ignore it.)
    
    return False, "No primary alcohol group found"


# Optional: Testing the function on a subset of example SMILES.
if __name__ == "__main__":
    # List a few example (Name, SMILES) pairs.
    test_examples = [
        ("(2S)-2-aminobutan-1-ol(1+)", "CC[C@@H](CO)[NH3+]"),
        ("Semiplenamide A", "O=C(NCCO)/C(=C/CC/C=C/CCCCCCCCCCCCC)/C"),
        ("3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,10-heptadecafluoro-1-decanol", "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
        ("furfuryl alcohol", "OCc1ccco1"),
        ("Ethanolamine", "NCCO"),
        ("docosan-1-ol", "OCCCCCCCCCCCCCCCCCCCCCO"),
        ("2-penten-1-ol", "OCC(=C(CC)[H])[H]")
    ]
    
    for name, smi in test_examples:
        result, reason = is_primary_alcohol(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")