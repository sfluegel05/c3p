"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: Anthoxanthin – flavonoid pigments in plants.
Anthoxanthins are typically based on a flavone (2-phenylchromen-4-one) scaffold 
but are often extensively substituted (e.g., with methoxy, hydroxy, and glycosyl groups)
that impart water solubility and the characteristic creamy to yellow color.
This implementation:
  • Excludes known false positives
  • Searches for a relaxed flavone (chromen-4-one) core using SMARTS.
  • Checks that the core carries an extra aromatic substituent (the B-ring of a flavone).
  • Ensures the molecular weight is in a reasonable range (240–1000 Da).
  • Requires the presence of oxygenated substituents (at least one methoxy or at least two hydroxy groups).
If any error occurs during parsing or matching, the function returns (None, None).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin (a plant flavonoid pigment) based on its SMILES string.
    
    The function:
      1. Parses the given SMILES.
      2. Applies an exclusion filter (e.g., for a known false positive such as Trichophenol A).
      3. Searches for a relaxed flavone core (chromen-4-one) using a revised SMARTS pattern.
      4. Verifies that the flavone core bears an attached aromatic substituent (the B-ring in 2-phenylchromen-4-one).
      5. Checks that molecular weight is within a typical range (240–1000 Da).
      6. Verifies that the molecule has enough oxygenated substituents 
         (at least one methoxy group or at least two hydroxy groups).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an anthoxanthin, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclusion filter: known false positive (for example, Trichophenol A)
    exclusion_smarts = "O=C1OC(C2=C(O)C=C(O)C=C2C)=CC=3C1=C(O)C=C(O)C3"
    exclusion_pattern = Chem.MolFromSmarts(exclusion_smarts)
    if exclusion_pattern and mol.HasSubstructMatch(exclusion_pattern):
        return False, "Matches exclusion pattern for Trichophenol A, not an anthoxanthin"
    
    # Look for the flavone core.
    # Revised SMARTS: "c1ccc2c(c1)oc(=O)c(c2)"
    # This pattern captures the benzopyran-4-one scaffold found in typical anthoxanthins.
    flavone_core_smarts = "c1ccc2c(c1)oc(=O)c(c2)"
    flavone_core = Chem.MolFromSmarts(flavone_core_smarts)
    if flavone_core is None:
        return None, None  # SMARTS parsing error
    
    core_matches = mol.GetSubstructMatches(flavone_core)
    if not core_matches:
        return False, "Flavone core (chromen-4-one) not found"
    
    # Check for the B-ring: a typical 2-phenyl substituent should be attached to the flavone core.
    # We search the neighbors of atoms in the matched core for an external aromatic carbon that belongs to a 6-membered ring.
    found_b_ring = False
    for match in core_matches:
        match_atoms = set(match)
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Skip atoms already part of the core
                if nbr_idx in match_atoms:
                    continue
                # Look for an aromatic carbon attached outside the core
                if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
                    # Check if this neighbor belongs to a 6-membered ring
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
    
    # Check molecular weight is within a typical range for anthoxanthins (240–1000 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 240:
        return False, f"Molecular weight {mol_wt:.1f} Da is too low for anthoxanthin"
    if mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} Da is too high for anthoxanthin"
    
    # Check for oxygenated substituents.
    # Count methoxy groups using SMARTS pattern "[OX2][CH3]"
    methoxy_pattern = Chem.MolFromSmarts("[OX2][CH3]")
    methoxy_count = len(mol.GetSubstructMatches(methoxy_pattern)) if methoxy_pattern else 0
    # Count hydroxy groups using SMARTS pattern "[OX2H]"
    hydroxy_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxy_count = len(mol.GetSubstructMatches(hydroxy_pattern)) if hydroxy_pattern else 0
    # Require at least one methoxy group or at least two hydroxy groups typical for anthoxanthins.
    if methoxy_count < 1 and hydroxy_count < 2:
        return False, "Insufficient oxygenated substituents (methoxy or hydroxy) typical of anthoxanthins"
    
    return True, "Molecule contains a flavone core with an attached B-ring and appropriate oxygenated substituents typical of anthoxanthins"

# Example usage (testing):
if __name__ == '__main__':
    # Test examples from the provided list
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