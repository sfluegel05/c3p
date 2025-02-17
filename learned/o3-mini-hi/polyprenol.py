"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: Polyprenol
Definition:
  A polyprenol is any prenol with the general formula H-[CH2C(Me)=CHCH2]nOH (n > 1).
  Hence, the molecule must have a terminal primary alcohol group (i.e. CH2OH where
  the CH2 is only connected to one heavy atom) and a carbon backbone composed of at
  least two consecutive isoprene units.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    
    A polyprenol must have:
      1) a terminal primary alcohol group (–CH2OH, where the CH2 is attached only to one heavy atom),
      2) at least one occurrence of two consecutive isoprene units, where an isoprene unit is defined as CH2–C(CH3)=CH–CH2.
         The detection is performed with the SMARTS pattern for two connected isoprene units:
         "[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]".
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple with True and a reason if the molecule is a polyprenol,
                     or False and the reason why it fails the classification.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # IMPORTANT: Add explicit hydrogens to ensure proper detection of terminal -CH2OH.
    molH = Chem.AddHs(mol)
    
    # STEP 1: Check for a terminal primary alcohol group.
    # Loop over all oxygen atoms in the hydrogen-added molecule.
    found_terminal_alcohol = False
    for atom in molH.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Oxygen
            # Check if this oxygen is part of an –OH (only one neighbor, which should be carbon).
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1:
                neighbor = neighbors[0]
                if neighbor.GetAtomicNum() == 6:  # Carbon neighbor
                    # Check whether this carbon is "primary" (only one heavy neighbor besides the oxygen)
                    heavy_neighbors = [nbr for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() != 1 and nbr.GetIdx() != atom.GetIdx()]
                    if len(heavy_neighbors) == 1:
                        found_terminal_alcohol = True
                        break
    if not found_terminal_alcohol:
        return False, "Molecule does not contain a terminal primary alcohol group"
    
    # STEP 2: Check for at least one block of two consecutive isoprene units.
    # The isoprene repeating unit is defined as CH2–C(CH3)=CH–CH2.
    # Two consecutive units yield the following SMARTS pattern:
    double_isoprene_pattern = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]")
    ip_matches = mol.GetSubstructMatches(double_isoprene_pattern)
    
    # Filter the matches so that every atom in the match is carbon.
    valid_matches = []
    for match in ip_matches:
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in match):
            valid_matches.append(match)
    
    if not valid_matches:
        return False, "No consecutive isoprene units found (at least 2 required)"
    
    # (Optional) Minimal molecular weight check.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} Da is too low for a polyprenol"
    
    # Passed both checks for terminal alcohol and isoprene units.
    n_blocks = len(valid_matches)
    return True, f"Contains terminal primary alcohol and {n_blocks} occurrence(s) of consecutive isoprene units."

# Example test calls (this section can be removed or commented out in production)
if __name__ == "__main__":
    test_examples = {
        "Glisoprenin F": "OC(CCC(O)C(O)(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)C)C)C)C)C",
        "ditrans,polycis-undecaprenol": "C(/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)/C)/C)/C)/C)/C)/C)/C)",
        "solanesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "(2E,6E,10E)-omega-hydroxyfarnesol": "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "(2-cis,6-trans)-farnesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CO",
        "geraniol": "CC(C)=CCC\\C(C)=C\\CO"
    }
    for name, sm in test_examples.items():
        res, reason = is_polyprenol(sm)
        print(f"Name: {name}\nSMILES: {sm}\n-> polyprenol: {res}\nReason: {reason}\n")