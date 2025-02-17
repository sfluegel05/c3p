"""
Classifies: CHEBI:15734 primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Primary Alcohol

Definition:
  A primary alcohol is a compound in which a hydroxyl group (-OH)
  is attached to a saturated (sp3) carbon atom. The carbon must be 
  primary – that is, it is either a methanol-type carbon that has 3 hydrogens,
  or a CH2OH-type where the carbon is attached to one other carbon (and 2 hydrogens).
  
Algorithm improvements over the previous version:
  - We add explicit hydrogens so that we can count them.
  - We inspect every oxygen atom and only consider those that (a) have at least one hydrogen,
    (b) are attached to exactly one heavy (non-H) atom.
  - For that heavy neighbor (which must be carbon and sp3),
    we count its heavy neighbors excluding the current –OH oxygen and count its hydrogens.
    If it is methanol-type (no other heavy neighbors and three H's) or CH2OH-type (one heavy neighbor,
    that neighbor must be carbon, and two H's) then we return True.
"""

from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule contains a primary alcohol group based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if a primary alcohol group is detected,
                     and False otherwise. The second element gives a reason for the classification.
    """
    # Convert SMILES to RDKit molecule and add explicit hydrogens so we have counts available.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Iterate over all atoms searching for candidate oxygens in –OH groups.
    for oxygen in mol.GetAtoms():
        if oxygen.GetSymbol() != "O":
            continue
        
        # Check if oxygen has at least one hydrogen attached.
        h_neighbors = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetSymbol() == "H"]
        if len(h_neighbors) < 1:
            continue  # Not an –OH group.
            
        # The oxygen should be attached to exactly one heavy (non-hydrogen) atom.
        heavy_neighbors = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) != 1:
            continue  # Avoid groups that are not classical alcohols (e.g. carboxyls or ethers).
        
        # The only heavy neighbor must be a carbon.
        carbon = heavy_neighbors[0]
        if carbon.GetSymbol() != "C":
            continue
        
        # The candidate carbon must be saturated (sp3 hybridized)
        if carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Count the carbon’s attached hydrogens using RDKit’s explicit hydrogens.
        # Note: After AddHs, all hydrogens are explicit.
        carbon_h_count = sum(1 for nbr in carbon.GetNeighbors() if nbr.GetSymbol() == "H")
        
        # Also count heavy neighbors attached to the candidate carbon (excluding the current –OH oxygen).
        heavy_nbrs_carbon = [nbr for nbr in carbon.GetNeighbors() 
                             if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != oxygen.GetIdx()]
        
        # Case A: methanol-type (CH3OH) where the candidate carbon has no other heavy neighbor
        # and exactly three hydrogens.
        if len(heavy_nbrs_carbon) == 0 and carbon_h_count == 3:
            return True, "Found primary alcohol group (methanol-type: CH3OH)"
        
        # Case B: CH2OH-type where the candidate carbon has exactly one heavy neighbor (and that neighbor must be C)
        # and exactly two hydrogens.
        if len(heavy_nbrs_carbon) == 1 and heavy_nbrs_carbon[0].GetSymbol() == "C" and carbon_h_count == 2:
            return True, "Found primary alcohol group (CH2OH-type)"
    
    # If we have not returned by now then no primary alcohol group was detected.
    return False, "No primary alcohol group found"

# Optional test suite (run as script)
if __name__ == "__main__":
    # Example test cases (name, SMILES) demonstrating positive examples.
    test_examples = [
        ("(2S)-2-aminobutan-1-ol(1+)", "CC[C@@H](CO)[NH3+]"),
        ("Semiplenamide A", "O=C(NCCO)/C(=C/CC/C=C/CCCCCCCCCCCCC)/C"),
        ("3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,10-heptadecafluoro-1-decanol",
         "OCCC(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
        ("(E,E)-2-methyl-6-oxohepta-2,4-dienol", "CC(=O)\\C=C\\C=C(/C)CO"),
        ("Pentadecanoyl-EA", "O=C(NCCO)CCCCCCCCCCCCCC"),
        ("(11Z)-icos-11-en-1-ol", "CCCCCCCC\\C=C/CCCCCCCCCCO"),
        ("syringin", "COc1cc(cc(OC)c1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)\\C=C\\CO"),
        ("2-(2-hydroxy-4-methoxyphenyl)-5-(3-hydroxypropyl)benzofuran",
         "COc1ccc(-c2cc3cc(CCCO)ccc3o2)c(O)c1"),
        ("2-penten-1-ol", "OCC(=C(CC)[H])[H]"),
        ("docosan-1-ol", "OCCCCCCCCCCCCCCCCCCCCCO"),
        ("ethanolamine", "NCCO")
    ]
    
    for name, smi in test_examples:
        result, reason = is_primary_alcohol(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")