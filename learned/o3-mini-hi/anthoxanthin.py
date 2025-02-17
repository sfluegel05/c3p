"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin – flavonoid pigments in plants.
Anthoxanthins are typically based on a flavone (2-phenylchromen-4-one) scaffold 
but are often extensively substituted (e.g., with methoxy, hydroxy, and glycosyl groups)
that impart water solubility and the characteristic creamy to yellow color.
This implementation:
  • Searches for a relaxed flavone (chromen-4-one) core using SMARTS.
  • Checks that the core carries an extra aromatic substituent (the B-ring of a flavone).
  • Ensures the molecular weight is in a reasonable range (240–1000 Da).
  • Requires the presence of oxygenated substituents (at least one methoxy or at least two hydroxy groups).
  • Excludes a known false positive (Trichophenol A) based on its pattern.
If any error occurs during parsing or matching, the function returns (None, None).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin (a plant flavonoid pigment) based on its SMILES string.
    
    The function:
      1. Parses the given SMILES.
      2. Applies an exclusion filter (e.g., for Trichophenol A).
      3. Searches for a relaxed flavone core (chromen-4-one) using SMARTS.
      4. Verifies that the flavone core bears an external aromatic substituent 
         (the phenyl B-ring required for a 2-phenylchromen-4-one scaffold).
      5. Checks that molecular weight is within a typical range (240–1000 Da).
      6. Verifies that the molecule has enough oxygenated substituents 
         (at least one methoxy group or at least two hydroxy groups).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an anthoxanthin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclusion filter: known false positive (Trichophenol A)
    exclusion_smarts = "O=C1OC(C2=C(O)C=C(O)C=C2C)=CC=3C1=C(O)C=C(O)C3"
    exclusion_pattern = Chem.MolFromSmarts(exclusion_smarts)
    if exclusion_pattern and mol.HasSubstructMatch(exclusion_pattern):
        return False, "Matches exclusion pattern for Trichophenol A, not an anthoxanthin"
    
    # Look for a flavone core. Here a relaxed SMARTS "c1cc2oc(=O)cc2c1" is used
    # which represents a chromen-4-one substructure (fused benzene and pyranone rings).
    flavone_core_smarts = "c1cc2oc(=O)cc2c1"
    flavone_core = Chem.MolFromSmarts(flavone_core_smarts)
    if flavone_core is None:
        return None, None  # parsing error for the core SMARTS
    
    core_matches = mol.GetSubstructMatches(flavone_core)
    if not core_matches:
        return False, "Flavone core (chromen-4-one) not found"

    # In a typical flavone (2-phenylchromen-4-one), a phenyl (B-ring) should be attached to the core.
    # We now check that at least one atom in a core match has a neighboring aromatic carbon that
    # is part of a 6-membered ring (a hint for a benzene ring).
    found_b_ring = False
    for match in core_matches:
        match_set = set(match)
        # Iterate over each atom in the core substructure match.
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Consider only neighbors not part of the flavone core that are aromatic carbons.
                if nbr_idx not in match_set and nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
                    # Check if this neighbor is in a 6-membered ring.
                    for ring in mol.GetRingInfo().AtomRings():
                        if nbr_idx in ring and len(ring) == 6:
                            found_b_ring = True
                            break
                if found_b_ring:
                    break
            if found_b_ring:
                break
        if found_b_ring:
            break
    if not found_b_ring:
        return False, "Flavone core found but no attached aromatic (B-ring) detected"
    
    # Molecular weight filter:
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 240:
        return False, f"Molecular weight {mol_wt:.1f} Da is too low for anthoxanthin"
    if mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} Da is too high for anthoxanthin"
    
    # Check for oxygenated substituents.
    # Count methoxy groups (-OCH3) using a SMARTS pattern.
    methoxy_pattern = Chem.MolFromSmarts("[OX2][CH3]")
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern)) if methoxy_pattern else 0
    # Count hydroxy groups (-OH) using a SMARTS pattern.
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern)) if hydroxy_pattern else 0
    # Require at least one methoxy or at least two hydroxy groups for typical anthoxanthins.
    if methoxy_count < 1 and hydroxy_count < 2:
        return False, "Insufficient oxygenated substituents (methoxy or hydroxy) typical of anthoxanthins"
    
    return True, "Molecule contains a flavone core with an attached B-ring and appropriate oxygenated substituents typical of anthoxanthins"

# Example usage (testing):
if __name__ == '__main__':
    # Some test examples (names and SMILES) from the list.
    test_smiles = {
        "sinensetin": "COc1ccc(cc1OC)-c1cc(=O)c2c(OC)c(OC)c(OC)cc2o1",
        "tambulin": "COc1ccc(cc1)-c1oc2c(OC)c(OC)cc(O)c2c(=O)c1O",
        "wogonin": "COc1c(O)cc(O)c2c1oc(cc2=O)-c1ccccc1",
        "xanthomicrol": "COc1c(O)c2c(oc(cc2=O)-c2ccc(O)cc2)c(OC)c1OC",
        "Trichophenol A": "O=C1OC(C2=C(O)C=C(O)C=C2C)=CC=3C1=C(O)C=C(O)C3",
        "6-geranylchrysin": "CC(C)=CCC\\C(C)=C\\Cc1c(O)cc2oc(cc(=O)c2c1O)-c1ccccc1"
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_anthoxanthin(smi)
        print(f"Name: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n")