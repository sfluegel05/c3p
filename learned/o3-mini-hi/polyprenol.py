"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: Polyprenol
Definition:
  A polyprenol is any prenol possessing the general formula:
       H-[CH2C(Me)=CHCH2]nOH with n>1,
  i.e. it must contain at least two isoprene units and a terminal primary alcohol
  (a –CH2OH group at one end of the carbon chain).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    
    Requirements:
      1) It must contain a terminal primary alcohol group.
         Here we search for an -OH group attached to a CH2 that is terminal. 
         We define "terminal" by requiring that the carbon is sp3 (tetrahedral),
         carries exactly two hydrogens, and has exactly one heavy-atom neighbor aside from the –OH.
      2) It must contain at least one occurrence of two consecutive isoprene units.
         An isoprene unit is defined as CH2–C(CH3)=CH–CH2.
         Two consecutive isoprene units are determined using the SMARTS:
           "[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]"
      3) Optionally, a reasonable molecular weight (>200 Da) is required.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple containing True and a descriptive reason if classified as a polyprenol,
                     or False with a reason explaining why it fails.
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for proper substructure detection.
    molH = Chem.AddHs(mol)
    
    # STEP 1: Check for a terminal primary alcohol group.
    # Instead of a strict SMARTS we iterate over oxygen atoms.
    terminal_alcohol_found = False
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 8:  # oxygen atom
            # Check that oxygen is part of an –OH group: one hydrogen neighbor.
            o_neighbors = atom.GetNeighbors()
            if sum(1 for nbr in o_neighbors if nbr.GetAtomicNum() == 1) != 1:
                continue
            # Get the heavy-atom neighbor (should be a carbon)
            carbon_neighbors = [nbr for nbr in o_neighbors if nbr.GetAtomicNum() == 6]
            if len(carbon_neighbors) != 1:
                continue
            c_atom = carbon_neighbors[0]
            # Check that this carbon is sp3 (tetrahedral) and has exactly two hydrogens.
            # Use explicit hydrogens (already added) to count.
            if c_atom.GetHybridization().name != "SP3":
                continue
            # Count hydrogen atoms attached to carbon.
            h_count = sum(1 for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() == 1)
            if h_count != 2:
                continue
            # Check terminality: carbon should have exactly one heavy neighbor aside from the -OH oxygen.
            heavy_neighbors = [nbr for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != atom.GetIdx()]
            if len(heavy_neighbors) == 1:
                terminal_alcohol_found = True
                break
    if not terminal_alcohol_found:
        return False, "Molecule does not contain a terminal primary alcohol group at a chain terminus"
    
    # STEP 2: Look for at least one occurrence of two consecutive isoprene units.
    # The SMARTS pattern is defined to detect CH2C(CH3)=CHCH2C(CH3)=CHCH2 segments.
    double_isoprene_pattern = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]")
    if double_isoprene_pattern is None:
        return False, "Error constructing consecutive isoprene units SMARTS pattern"
    
    isoprene_matches = mol.GetSubstructMatches(double_isoprene_pattern)
    if not isoprene_matches:
        return False, "No consecutive isoprene units found (min. two consecutive isoprene units required)"
    
    # OPTIONAL: Check that the molecule has a reasonable molecular weight (>200 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a polyprenol"
    
    n_occurrences = len(isoprene_matches)
    return True, f"Contains terminal primary alcohol and {n_occurrences} occurrence(s) of consecutive isoprene units."

# Example test calls (these may be commented out in production)
if __name__ == "__main__":
    test_examples = {
      "Glisoprenin F": "OC(CCC(O)C(O)(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)C)C)C)C)C",
      "ditrans,polycis-undecaprenol": "C(/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)/C)/C)/C)/C)/C)/C)/C)",
      "solanesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
      "(2E,6E,10E)-omega-hydroxyfarnesol": "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",
      "(2-cis,6-trans)-farnesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CO",
      "all-trans-hexaprenol": "C(/C=C(/CCC=C(C)C)\\C)C\\C(\\C)=C\\CC\\C(\\C)=C\\CC/C(/C)=C/CC\\C(=C\\CO)\\C",
      "geraniol": "CC(C)=CCC\\C(C)=C\\CO",
      # Further examples can be added similarly.
    }
    for name, sm in test_examples.items():
        res, reason = is_polyprenol(sm)
        print(f"Name: {name}\nSMILES: {sm}\n-> polyprenol: {res}\nReason: {reason}\n")