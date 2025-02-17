"""
Classifies: CHEBI:26199 polyprenol
"""
"""
Classifies: Polyprenol [Any member of the class of prenols possessing the general formula 
H-[CH2C(Me)=CHCH2]nOH, with n > 1 (more than one isoprene unit)].
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    
    A polyprenol must have a terminal primary alcohol group 
    (a CH2–OH that is attached to only one heavy atom) and 
    a carbon backbone composed of at least two consecutive isoprene units.
    The isoprene repeating unit is represented as CH2–C(CH3)=CH–CH2.
    We enforce contiguity by requiring a match to a pattern with two
    consecutive isoprene units, i.e. "[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]".
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple with the classification result and the reason.
    """
    
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # STEP 1: Check for a terminal primary alcohol.
    # First get all matches for a primary alcohol: a CH2 linked with –OH.
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2][OX2H]")
    pa_matches = mol.GetSubstructMatches(primary_alcohol_pattern)
    has_terminal_alcohol = False
    for match in pa_matches:
        # match is a tuple of atom indices; the first atom in our pattern is the CH2.
        ch2_idx = match[0]
        ch2_atom = mol.GetAtomWithIdx(ch2_idx)
        # Count heavy-atom neighbors (neighbors that are not hydrogen)
        heavy_neighbors = [nbr for nbr in ch2_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) == 1:
            has_terminal_alcohol = True
            break
    if not has_terminal_alcohol:
        return False, "Molecule does not contain a terminal primary alcohol group"
    
    # STEP 2: Check for at least two consecutive isoprene units.
    # Define a SMARTS pattern for a block of two consecutive isoprene units.
    # Isoprene repeating unit: CH2–C(CH3)=CH–CH2.
    # Two consecutive units: [CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]
    double_isoprene_pattern = Chem.MolFromSmarts("[CH2]C([CH3])=C[CH2]C([CH3])=C[CH2]")
    ip_matches = mol.GetSubstructMatches(double_isoprene_pattern)
    
    # Filter matches so that every atom in the match is carbon.
    valid_matches = []
    for match in ip_matches:
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in match):
            valid_matches.append(match)
    
    if not valid_matches:
        return False, "No consecutive isoprene units found (at least 2 required)"
    
    # Optionally, check a minimal molecular weight cutoff (e.g. 200 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight {mol_wt:.1f} Da is too low for a polyprenol"
    
    # If both conditions are met, then we classify as a polyprenol.
    # We report how many consecutive isoprene blocks were found.
    n_blocks = len(valid_matches)
    return True, f"Contains terminal primary alcohol and {n_blocks} occurrence(s) of two consecutive isoprene units."

# Example test calls (this part may be removed in production if only the function is needed)
if __name__ == "__main__":
    test_examples = {
        "Glisoprenin F": "OC(CCC(O)C(O)(C)C)(CCCC(O)(CCCC(O)(CCCC(O)(CCCC(O)(CCC/C(=C/CC/C(=C/CC/C(=C/CO)/C)/C)/C)C)C)C)C)C",
        "ditrans,polycis-undecaprenol": "C(/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(\\CC/C=C(/CC/C=C(/CCC=C(C)C)\\C)\\C)/C)/C)/C)/C)/C)/C)/C)",
        "solanesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "(2E,6E,10E)-omega-hydroxyfarnesol": "C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CO",
        "(2-cis,6-trans)-farnesol": "CC(C)=CCC\\C(C)=C\\CC\\C(C)=C/CO",
        "geraniol (prenol)": "CC(C)=CCC\\C(C)=C\\CO"
    }
    for name, sm in test_examples.items():
        res, reason = is_polyprenol(sm)
        print(f"Name: {name}\nSMILES: {sm}\n-> polyprenol: {res}\nReason: {reason}\n")