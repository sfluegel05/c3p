"""
Classifies: CHEBI:18254 ribonucleoside
"""
"""
Classifies: any nucleoside where the sugar component is D-ribose (ribonucleoside)
Definition: A ribonucleoside has a nucleobase linked to a D-ribofuranose unit.
Note: This implementation uses two approaches. First, it checks for common SMARTS
patterns for an unmodified sugar ([C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O) or for a 2'-O-methyl modified sugar (CO[C@@H]1[C@H](O)[C@@H](CO)O[C@H]1).
If these do not match, the code then scans for any five-membered ring with exactly one oxygen and four carbons
and further looks for an exocyclic -CH2OH substituent on one of the carbons – a hallmark of D-ribose.
Because nucleosides are highly variable, this function is heuristic in nature.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ribonucleoside(smiles: str):
    """
    Determines if a molecule is a ribonucleoside (a nucleoside with a D-ribose sugar)
    based on its SMILES string.
    
    The function uses two strategies:
      1. It checks for common SMARTS fragments found in ribonucleosides:
         - an unmodified ribose pattern: "[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"
         - a 2'-O-methyl modified sugar pattern: "CO[C@@H]1[C@H](O)[C@@H](CO)O[C@H]1"
      2. If no match is found, it searches the molecule for any five-membered ring that
         contains exactly one oxygen and four carbons and then checks if a ring carbon bears
         an exocyclic CH2OH fragment.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a ribonucleoside, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First approach: attempt to match common SMARTS patterns for a D-ribose sugar.
    # Pattern 1: Unmodified D-ribofuranose unit (as seen in many canonical ribonucleosides)
    ribose_smarts1 = "[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O"
    pattern1 = Chem.MolFromSmarts(ribose_smarts1)
    if pattern1 is None:
        return False, "Error creating SMARTS pattern for ribose (pattern1)"
    
    # Pattern 2: A variant (e.g., for 2'-O-methyl nucleosides)
    ribose_smarts2 = "CO[C@@H]1[C@H](O)[C@@H](CO)O[C@H]1"
    pattern2 = Chem.MolFromSmarts(ribose_smarts2)
    if pattern2 is None:
        return False, "Error creating SMARTS pattern for ribose (pattern2)"
    
    if mol.HasSubstructMatch(pattern1):
        return True, "Matched common D-ribose pattern (unmodified ribofuranose)"
    if mol.HasSubstructMatch(pattern2):
        return True, "Matched common D-ribose pattern with 2'-O-methyl modification"
    
    # Second approach: try a heuristic analysis of the rings.
    # Add hydrogens for better analysis of substituents.
    mol_h = Chem.AddHs(mol)
    ring_info = mol_h.GetRingInfo()
    # Iterate over all rings in the molecule
    for ring in ring_info.AtomRings():
        if len(ring) == 5:
            # Count how many atoms in the ring are oxygen and how many are carbons
            oxy_count = 0
            carbon_count = 0
            for idx in ring:
                atom = mol_h.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 8:
                    oxy_count += 1
                elif atom.GetAtomicNum() == 6:
                    carbon_count += 1
            # For a ribofuranose ring, we expect exactly one oxygen and four carbons.
            if oxy_count == 1 and carbon_count == 4:
                # Look for exocyclic CH2OH group on one of the ring carbons.
                for idx in ring:
                    atom = mol_h.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() == 6:  # Consider only carbon atoms
                        for nbr in atom.GetNeighbors():
                            if nbr.GetIdx() not in ring and nbr.GetAtomicNum() == 6:
                                # Check if the substituent carbon (nbr) could be CH2OH:
                                # It should be attached to an oxygen (that is not in the ring)
                                oxy_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() == 8 and n.GetIdx() not in ring]
                                # Also, require that the substituent carbon has at least two hydrogens.
                                if oxy_neighbors and nbr.GetTotalNumHs() >= 2:
                                    return True, "Found five‐membered ring with exocyclic CH2OH substituent consistent with D‐ribose"
    return False, "No D-ribose moiety (ribonucleoside sugar) found"