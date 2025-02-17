"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: Inositol Phosphoceramide
Definition: A phosphosphingolipid in which an inositol residue and the ceramide moiety are 
linked via a phosphodiester bridge. The ceramide moiety contains substituents (typically an amide bond
and long aliphatic chains) that vary with different sphingoid bases and fatty acyl moieties.
This implementation:
  - Searches for an inositol ring with a relaxed SMARTS (allowing alternate stereo annotations).
  - Looks for an amide bond (C(=O)N) as evidence of the ceramide part.
  - Iterates phosphorus atoms to find one that has at least one P=O and at least two single-bonded oxygens,
    with the single bonds connecting one branch to the inositol ring and one branch to non-inositol (the ceramide).
  - Requires the molecule to meet minimal molecular weight and carbon count.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    The method checks:
      1) The presence of an inositol ring substructure: a cyclohexane bearing an OH on each carbon.
         (Here we try a main SMARTS and an alternate as backup.)
      2) The presence of a ceramide moiety indicated by an amide bond (C(=O)N).
      3) A phosphodiester linkage in which a phosphorus atom (P) has at least one double‐bonded oxygen (P=O)
         and two single‐bonded oxygens, among which one branch is linked to the inositol substructure
         and the other branch goes to the ceramide side.
      4) A minimum overall molecular weight and number of carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an inositol phosphoceramide, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Find the inositol ring ---
    # Try a primary SMARTS pattern for myo-inositol as seen in many compounds.
    inositol_smarts1 = "O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"
    inositol_query = Chem.MolFromSmarts(inositol_smarts1)
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    # If no match, try an alternate pattern that is less stereochemically strict.
    if not inositol_matches:
        inositol_smarts2 = "OC1C(O)C(O)C(O)C(O)C1O"
        inositol_query = Chem.MolFromSmarts(inositol_smarts2)
        inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        return False, "Missing inositol ring substructure"
    
    # Collect all atom indices that are part of an inositol match.
    inositol_atom_indices = set()
    for match in inositol_matches:
        inositol_atom_indices.update(match)
    
    # --- Step 2: Find evidence of a ceramide moiety (look for at least one amide bond: C(=O)N) ---
    amide_smarts = "C(=O)N"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    amide_matches = mol.GetSubstructMatches(amide_query)
    if not amide_matches:
        return False, "Missing ceramide amide bond"
    
    # For reference we collect the carbonyl carbons involved in amide bonds.
    amide_carbons = set(match[0] for match in amide_matches)
    
    # --- Step 3: Check for the phosphodiester linkage ---
    phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:   # phosphorus atomic number
            continue
        p_atom = atom
        bonds = p_atom.GetBonds()
        # Collect oxygens attached to phosphorus via different bond types.
        double_bonded_os = []
        single_bonded_os = []
        for bond in bonds:
            nbr = bond.GetOtherAtom(p_atom)
            if nbr.GetAtomicNum() != 8:
                continue
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                double_bonded_os.append(nbr)
            elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                single_bonded_os.append(nbr)
        # Require at least one P=O and at least two single-bonded oxygens.
        if len(double_bonded_os) < 1 or len(single_bonded_os) < 2:
            continue
        
        # Check the connectivity: 
        # - Look for at least one single-bonded oxygen that connects (via its other bond) 
        #   to the inositol ring.
        # - Look for at least one single-bonded oxygen that connects to a fragment outside inositol.
        inositol_connected = False
        ceramide_connected = False
        for oxy in single_bonded_os:
            # For each oxygen, check its neighbors (other than phosphorus)
            for nbr in oxy.GetNeighbors():
                if nbr.GetIdx() == p_atom.GetIdx():
                    continue
                # If the neighbor is part of the inositol substructure, mark it.
                if nbr.GetIdx() in inositol_atom_indices:
                    inositol_connected = True
                else:
                    # As a proxy for ceramide connectivity, we require either
                    # that atom is not in the inositol set OR it is near an amide bond.
                    ceramide_connected = True
            if inositol_connected and ceramide_connected:
                phosphate_found = True
                break
        if phosphate_found:
            break
    
    if not phosphate_found:
        return False, "Phosphodiester linkage (connecting phosphate to both inositol and ceramide) not detected"
    
    # --- Step 4: Overall molecular size check ---
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for inositol phosphoceramide"
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Too few carbons ({carbon_count}) to be an inositol phosphoceramide"
    
    return True, "Molecule contains an inositol ring, an amide bond (indicative of ceramide), and a phosphate group that bridges inositol to ceramide"

# Example usage for testing:
if __name__ == "__main__":
    # A sample SMILES from one of the provided examples:
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)CCCCCCCCCCCCCC"
    result, reason = is_inositol_phosphoceramide(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)