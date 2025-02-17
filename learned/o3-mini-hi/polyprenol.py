"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: Polyprenol
Definition:
  A polyprenol is any prenol with the general formula H-[CH2C(Me)=CHCH2]nOH (n > 1).
  Hence, the molecule must have a terminal primary alcohol group (-CH2OH, where the CH2 is
  connected to exactly one non-hydrogen atom other than the hydroxyl oxygen), and it must 
  contain at least one block of two consecutive isoprene units.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    
    A polyprenol must have:
      1) A terminal primary alcohol group. This is detected by looking for a CH2 group
         bonded to an -OH such that the CH2 carries exactly two hydrogens (and therefore,
         only one heavy neighbor apart from the oxygen) and the oxygen is only bonded to that 
         CH2 (plus one hydrogen).
      2) At least one occurrence of two consecutive isoprene units.
         Each isoprene unit is defined as CH2–C(CH3)=CH–CH2.
         Two consecutive units match the SMARTS pattern: "[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]".
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple containing True and a reason if the molecule is a polyprenol,
                     or False with the reason why it fails the classification.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to ensure correct count for detecting the -CH2OH group.
    molH = Chem.AddHs(mol)
    
    # STEP 1: Check for a terminal primary alcohol group.
    # We iterate over carbon atoms and look for a CH2 that is bonded to an -OH.
    found_terminal_alcohol = False
    for atom in molH.GetAtoms():
        # Look for carbon atoms.
        if atom.GetAtomicNum() == 6:
            # Count total hydrogens attached (including implicit as explicit now).
            total_h = atom.GetTotalNumHs()
            # Gather heavy neighbors (atomic number > 1).
            heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
            # For a CH2 in a terminal primary alcohol, we expect:
            # - Exactly two hydrogens (CH2)
            # - Exactly two heavy neighbors: one being the oxygen of the OH, and the other being the chain.
            if total_h == 2 and len(heavy_neighbors) == 2:
                # Check if one of these heavy neighbors is oxygen and qualifies as -OH.
                for nbr in heavy_neighbors:
                    if nbr.GetAtomicNum() == 8:
                        # Check that the oxygen has exactly one heavy neighbor (the carbon we are examining)
                        o_heavy = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() > 1]
                        if len(o_heavy) == 1:
                            found_terminal_alcohol = True
                            break
            if found_terminal_alcohol:
                break
    
    if not found_terminal_alcohol:
        return False, "Molecule does not contain a terminal primary alcohol group"
    
    # STEP 2: Look for at least one occurrence of two consecutive isoprene units.
    # The SMARTS pattern for two consecutive isoprene units:
    double_isoprene_pattern = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]")
    if double_isoprene_pattern is None:
        return False, "Error in isoprene SMARTS pattern construction"
    # Find substructure matches in the original molecule (hydrogens not strictly needed here).
    isoprene_matches = mol.GetSubstructMatches(double_isoprene_pattern)
    if not isoprene_matches:
        return False, "No consecutive isoprene units found (at least 2 required)"
    
    # OPTIONAL: Check a minimal molecular weight to filter out fragments.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} Da is too low for a polyprenol"
    
    # Count the number of consecutive isoprene occurrences found.
    n_occurrences = len(isoprene_matches)
    
    return True, f"Contains terminal primary alcohol and {n_occurrences} occurrence(s) of consecutive isoprene units."

# Example test calls (can be removed or commented out in production)
if __name__ == "__main__":
    test_examples = {
        "Glisoprenin F": "OC(CCC(O)C(O)(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)C)C)C)C)C",
        "ditrans,polycis-undecaprenol": "C(/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)/C)/C)/C)/C)/C)/C)/C)",
        "solanesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "(2E,6E,10E)-omega-hydroxyfarnesol": "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "(2-cis,6-trans)-farnesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CO",
        "all-trans-hexaprenol": "C(/C=C(/CCC=C(C)C)\\C)C\\C(\\C)=C\\CC\\C(\\C)=C\\CC/C(/C)=C/CC\\C(=C\\CO)\\C",
        "geraniol": "CC(C)=CCC\\C(C)=C\\CO",
        # Additional examples can be added here.
    }
    for name, sm in test_examples.items():
        res, reason = is_polyprenol(sm)
        print(f"Name: {name}\nSMILES: {sm}\n-> polyprenol: {res}\nReason: {reason}\n")