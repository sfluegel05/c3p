"""
Classifies: CHEBI:15734 primary alcohol
"""
#!/usr/bin/env python3
"""
Classifies: Primary Alcohol

Definition:
  A primary alcohol is a compound in which a hydroxyl group (-OH)
  is attached to a saturated (sp3) carbon atom that is primary – that is,
  the carbon bears either:
    • three hydrogen atoms (as in methanol, CH3OH), or
    • two hydrogen atoms if it is attached to one (and only one) other carbon (CH2OH).

Algorithm improvements over the previous version:
  - We more strictly check that the –OH oxygen is attached to exactly one heavy atom.
  - For the candidate carbon we check its sp3 hybridization.
  - We count both the explicit hydrogen neighbors and the heavy neighbors (excluding the –OH oxygen).
    For a primary alcohol, the candidate carbon must have:
      • 0 other heavy neighbors and 3 hydrogens (methanol-type), or
      • 1 heavy neighbor (which must be a carbon) and 2 hydrogens (CH2OH-type).
  - This helps filter out cases where additional heavy atoms (e.g. from esters, carboxyls or secondary/tertiary alcohols)
    would otherwise be misinterpreted.
  
Note:
  RDKit’s AddHs method is used so that –OH hydrogens are available for counting.
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
    # Convert SMILES to RDKit molecule and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)

    # Loop over atoms to find candidate oxygen atoms.
    for atom in mol.GetAtoms():
        # We are interested in oxygen atoms.
        if atom.GetSymbol() != "O":
            continue
        
        # Check at least one explicit hydrogen attached to the oxygen.
        attached_h = [nbr for nbr in atom.GetNeighbors() if nbr.GetSymbol() == "H"]
        if len(attached_h) < 1:
            continue  # Not an –OH group if no hydrogen is attached.
        
        # The oxygen atom in an alcohol should be attached to exactly one heavy atom.
        heavy_nbrs = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_nbrs) != 1:
            continue  # Exclude cases like carboxylate or bridging ethers.
        
        # The only heavy neighbor should be a carbon.
        cand_carbon = heavy_nbrs[0]
        if cand_carbon.GetSymbol() != "C":
            continue
        
        # Candidate carbon must be sp3-hybridized (saturated).
        if cand_carbon.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Count explicit hydrogen atoms attached to the candidate carbon.
        carbon_h_count = sum(1 for nbr in cand_carbon.GetNeighbors() if nbr.GetSymbol() == "H")
        
        # Count heavy neighbors of the candidate carbon excluding the current –OH oxygen.
        heavy_nbrs_c = [nbr for nbr in cand_carbon.GetNeighbors() 
                        if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != atom.GetIdx()]
        
        # For a primary alcohol, the candidate carbon must be:
        #  a) methanol-type: no other heavy neighbor and exactly 3 hydrogens, or
        #  b) CH2OH-type: exactly one heavy neighbor (and that neighbor must be a carbon) and exactly 2 hydrogens.
        if len(heavy_nbrs_c) == 0 and carbon_h_count == 3:
            return True, "Found primary alcohol group (methanol-type: CH3OH)"
        elif len(heavy_nbrs_c) == 1 and heavy_nbrs_c[0].GetSymbol() == "C" and carbon_h_count == 2:
            return True, "Found primary alcohol group (CH2OH-type)"
        # Else, not a primary alcohol candidate based on this –OH.
    
    return False, "No primary alcohol group found"


# Below is an optional test suite – you can run this module to see example outputs.
if __name__ == "__main__":
    # List some test examples (Name, SMILES)
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