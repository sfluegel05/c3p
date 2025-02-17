"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: Polyprenol
Definition:
  A polyprenol is any prenol with the general formula H-[CH2C(Me)=CHCH2]nOH (n > 1).
  Thus, the molecule must contain a terminal primary alcohol group – that is, a CH2OH
  where the CH2 is connected to exactly one carbon (the chain) and the –OH – and
  at least one occurrence of two consecutive isoprene units.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    
    Requirements:
      1) Contains a terminal primary alcohol group.
         This is identified by the SMARTS pattern "[C;D2;H2][O;H1]". Here,
         [C;D2;H2] is a saturated (sp3) carbon atom having exactly two heavy-atom neighbors
         (degree 2) and exactly two hydrogens, and [O;H1] is an oxygen attached to one hydrogen.
         In addition, we check that the CH2’s non-oxygen neighbor is a terminal attachment 
         (i.e. it is the only heavy neighbor aside from the –OH oxygen).
      2) Contains at least one occurrence of two consecutive isoprene units.
         An isoprene unit is defined as CH2–C(CH3)=CH–CH2, and two consecutive
         isoprene units are detected using the SMARTS pattern:
           "[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]"
      3) (Optionally) The molecule should have a reasonable molecular weight (e.g. >200 Da)
          to help filter out fragments.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple containing True and a descriptive reason if the molecule is classified as 
                     a polyprenol, or False with a reason explaining why it fails.
    """
    # Parse molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for accurate substructure detection.
    molH = Chem.AddHs(mol)
    
    # STEP 1: Check for a terminal primary alcohol group using a SMARTS query.
    # The SMARTS "[C;D2;H2][O;H1]" matches a CH2 group (degree 2, 2 hydrogens) attached to an –OH group.
    alcohol_pattern = Chem.MolFromSmarts("[C;D2;H2][O;H1]")
    if alcohol_pattern is None:
        return False, "Error constructing terminal alcohol SMARTS pattern"
    
    matches = molH.GetSubstructMatches(alcohol_pattern)
    if not matches:
        return False, "Molecule does not contain a terminal primary alcohol group"
    
    # Ensure that at least one match is truly at the terminus of the chain.
    # We flag a match as terminal if the CH2 (first atom in the match) has a non-O neighbor count of exactly one.
    found_terminal = False
    for match in matches:
        # match returns indices: (carbon_index, oxygen_index)
        c_idx, o_idx = match
        c_atom = molH.GetAtomWithIdx(c_idx)
        # Get heavy neighbors other than the alcohol oxygen.
        non_O_neighbors = [nbr for nbr in c_atom.GetNeighbors() 
                           if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != o_idx]
        if len(non_O_neighbors) == 1:
            found_terminal = True
            break
    if not found_terminal:
        return False, "Molecule does not contain a terminal primary alcohol group at a chain terminus"
    
    # STEP 2: Look for at least one occurrence of two consecutive isoprene units.
    # SMARTS for two consecutive isoprene units:
    double_isoprene_pattern = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]")
    if double_isoprene_pattern is None:
        return False, "Error constructing consecutive isoprene units SMARTS pattern"
    
    isoprene_matches = mol.GetSubstructMatches(double_isoprene_pattern)
    if not isoprene_matches:
        return False, "No consecutive isoprene units found (at least 2 required)"
    
    # OPTIONAL: Filter on molecular weight to avoid small fragments.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} Da is too low for a polyprenol"
    
    n_occurrences = len(isoprene_matches)
    return True, f"Contains terminal primary alcohol and {n_occurrences} occurrence(s) of consecutive isoprene units."

# Example test calls (these may be commented out in production)
if __name__ == "__main__":
    test_examples = {
      "Glisoprenin F": "OC(CCC(O)C(O)(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)C)C)C)C)C",
      "ditrans,polycis-undecaprenol": "C(/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)/C)/C)/C)/C)/C)/C)/C)",
      "solanesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
      "(2E,6E,10E)-omega-hydroxyfarnesol": "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",
      "(2-cis,6-trans)-farnesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CO",
      "all-trans-hexaprenol": "C(/C=C(/CCC=C(C)C)\\C)C\\C(\\C)=C\\CC\\C(\\C)=C\\CC/C(/C)=C/CC\\C(=C\\CO)\\C",
      "geraniol": "CC(C)=CCC\\C(C)=C\\CO",
      # Additional examples omitted for brevity.
    }
    for name, sm in test_examples.items():
        res, reason = is_polyprenol(sm)
        print(f"Name: {name}\nSMILES: {sm}\n-> polyprenol: {res}\nReason: {reason}\n")