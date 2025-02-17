"""
Classifies: CHEBI:15734 primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Primary alcohol
Definition: A primary alcohol is a compound in which a hydroxyl group (-OH) is attached 
           to a saturated (sp3) carbon that is primary – having either no or only one other carbon attached.
           
In other words, the -OH must be bonded to a carbon that is either:
 • methanol-type: CH3–OH (the carbon has no other carbon neighbors)
 • CH2OH-type: R–CH2–OH (the carbon has exactly one carbon neighbor)
           
This improved algorithm first adds explicit hydrogens so that the true hydrogen count is visible.
It then checks every oxygen atom that bears a hydrogen, ensuring that it is attached to exactly one carbon.
Finally, the candidate carbon is verified to be sp³ and “primary” (by counting how many carbon neighbors it has excluding the OH oxygen).
"""

from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule contains a primary alcohol group based on its SMILES string.
    
    The algorithm works as follows:
      1. Parse the SMILES string and add explicit hydrogens.
      2. Loop over all oxygen atoms. For each oxygen that has at least one hydrogen,
         check that (a) it is attached to exactly one carbon atom (ensuring it is an -OH rather than a bridging oxygen),
         and (b) that carbon is sp³ and is primary (i.e. apart from the OH oxygen, it is attached to at most one other carbon).
      3. If any such –OH is found, return True along with a reason indicating whether it is a methanol-type or CH2OH-type.
      4. Otherwise, return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple of (True, reason) if a primary alcohol group is detected,
                     otherwise (False, explanation).
    """
    # Parse SMILES and add explicit hydrogens to reveal implicit H's.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Loop through all atoms to find candidate –OH groups.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "O":
            continue
        
        # Must have at least one hydrogen attached.
        if atom.GetTotalNumHs() < 1:
            continue

        # Get the heavy neighbors (ignoring hydrogen atoms)
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # For an –OH group we expect exactly one heavy (non‐H) neighbor.
        if len(heavy_neighbors) != 1:
            continue
        cand = heavy_neighbors[0]
        if cand.GetSymbol() != "C":
            continue

        # Make sure candidate carbon is sp3 hybridized.
        if cand.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Count how many carbon neighbors the candidate carbon has other than this OH oxygen.
        cand_carbon_neighbors = [nbr for nbr in cand.GetNeighbors() if nbr.GetSymbol()=="C" and nbr.GetIdx() != atom.GetIdx()]
        n_carbon_neighbors = len(cand_carbon_neighbors)
        # For a primary carbon (the one bearing the -OH), it should have at most one other carbon.
        if n_carbon_neighbors > 1:
            continue  # Not primary
        
        # (Optional) Verify hydrogen count on the candidate.
        # For CH3–OH (methanol-type) we expect 3 total H's; for R–CH2–OH we expect 2.
        h_count = cand.GetTotalNumHs()
        if n_carbon_neighbors == 0 and h_count == 3:
            return True, "Found primary alcohol group (methanol-type: CH3OH)"
        elif n_carbon_neighbors == 1 and h_count == 2:
            return True, "Found primary alcohol group (CH2OH-type)"
        # If the hydrogen count does not match exactly, we still accept the group (this may occur due to resonance or tautomers).
        elif n_carbon_neighbors <= 1:
            return True, "Found primary alcohol group (primary C–OH detected)"
    
    return False, "No primary alcohol group found"

# Optional: Some testing examples (remove or modify as needed)
if __name__ == "__main__":
    test_smiles = [
        "CC[C@@H](CO)[NH3+]",                           # (2S)-2-aminobutan-1-ol (CH2OH-type)
        "O=C(NCCO)/C(=C/CC/C=C/CCCCCCCCCCCCC)/C",        # Semiplenamide A (CH2OH-type)
        "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F",  # 1-decanol derivative (CH2OH-type)
        "OCc1ccco1",                                   # Furfuryl alcohol (CH2OH-type)
        "NCCO",                                        # Ethanolamine (CH2OH-type)
        "OCCCCCCCCCCCCCCCCCCCCCO",                     # docosan-1-ol (CH2OH-type)
        "OC1(CCC(CC1)C(C)=C)",                         # beta-terpineol (should be rejected as its OH is on a non‐primary center)
        "OC1C=C(Cl)C(Cl)C=C1Cl",                        # 2,4,5-trichlorocyclohexa-2,5-dien-1-ol (should be rejected)
        "C1=CC=C2C(=C1)[C@]34CCN5CC6=CCO[C@H]([C@]([C@]6(C[C@@]35[H])[H])([C@@]4(N2)[H])[H])O"  # caracurine VII (should be missed)
    ]
    
    for smi in test_smiles:
        result, reason = is_primary_alcohol(smi)
        print(f"SMILES: {smi}\nResult: {result}, Reason: {reason}\n")