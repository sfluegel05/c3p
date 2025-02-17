"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: Inositol Phosphoceramide
Definition: A phosphosphingolipid in which an inositol residue and the ceramide moiety are 
linked via a phosphodiester bridge. The ceramide moiety contains substituents (typically an amide bond
and long aliphatic chains) that vary with different sphingoid bases and fatty acyl moieties.
This implementation improves upon a previous attempt by:
  - Using two alternate SMARTS patterns to detect an inositol ring.
  - Requiring the presence of an amide bond (C(=O)N) as evidence of the ceramide.
  - Checking for a phosphorus atom bearing at least one P=O (double‐bonded oxygen) and two single‐bonded oxygens,
    one of which must be directly attached to an inositol atom and at least one branch (checked via a short bond path)
    must lead to an amide carbon.
  - Requiring the overall molecule to meet minimal molecular weight and carbon count thresholds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    The method checks:
      1) The presence of an inositol ring substructure using two alternate SMARTS patterns.
      2) The presence of a ceramide moiety indicated by an amide bond (C(=O)N).
      3) For at least one phosphorus atom (P):
            - It must have at least one P=O (double‐bonded oxygen) and at least two oxygens bound via single bonds.
            - Among its single‐bonded oxygens, one must attach to a carbon known to be part of an inositol ring.
            - At least one of the remaining oxygens must connect (via a short path, ≤3 bonds) to a carbonyl
              carbon from an amide bond.
      4) A minimum overall molecular weight and number of carbon atoms (to avoid very small molecules).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an inositol phosphoceramide, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- Step 1: Detect an inositol ring ---
    # First try a stereochemically strict myo-inositol pattern.
    inositol_smarts1 = "O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"
    inositol_query = Chem.MolFromSmarts(inositol_smarts1)
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        # If not found, try a less stereochemically strict pattern.
        inositol_smarts2 = "OC1C(O)C(O)C(O)C(O)C1O"
        inositol_query = Chem.MolFromSmarts(inositol_smarts2)
        inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        return False, "Missing inositol ring substructure"
    
    # Collect all atom indices that constitute the inositol(s)
    inositol_atom_indices = set()
    for match in inositol_matches:
        inositol_atom_indices.update(match)
    
    # --- Step 2: Find evidence for a ceramide moiety (look for an amide bond: C(=O)N) ---
    amide_smarts = "C(=O)N"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    amide_matches = mol.GetSubstructMatches(amide_query)
    if not amide_matches:
        return False, "Missing ceramide amide bond"
    
    # For later use, note the carbonyl carbons (first atom in the match) that appear in amide bonds.
    amide_carbons = set(match[0] for match in amide_matches)
    
    # --- Step 3: Check for a phosphodiester linkage that bridges inositol and ceramide ---
    phosphate_found = False
    # Loop over all phosphorus atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        p_atom = atom
        bonds = p_atom.GetBonds()
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
        
        found_inositol_branch = False
        found_ceramide_branch = False
        # For each single-bonded oxygen check connectivity.
        for oxy in single_bonded_os:
            # Look at neighbors of the oxygen besides the phosphorus.
            for nbr in oxy.GetNeighbors():
                if nbr.GetIdx() == p_atom.GetIdx():
                    continue
                # If the neighbor is part of the inositol ring, mark that branch.
                if nbr.GetIdx() in inositol_atom_indices:
                    found_inositol_branch = True
                else:
                    # Otherwise, check if this branch leads to a ceramide signature.
                    # We check if the current oxygen is within a short path (<=3 bonds) to any amide carbon.
                    for amide_c in amide_carbons:
                        path = Chem.rdmolops.GetShortestPath(mol, oxy.GetIdx(), amide_c)
                        if 0 < len(path) - 1 <= 3:
                            found_ceramide_branch = True
                            break
            # If already found both connections, no need to check further.
            if found_inositol_branch and found_ceramide_branch:
                phosphate_found = True
                break
        if phosphate_found:
            break

    if not phosphate_found:
        return False, "Phosphodiester linkage (bridging inositol and ceramide) not detected"
    
    # --- Step 4: Overall size checks ---
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for inositol phosphoceramide"
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Too few carbons ({carbon_count}) to be an inositol phosphoceramide"
    
    return True, "Molecule contains an inositol ring, an amide bond (ceramide), and a phosphate that bridges them"

# For example usage:
if __name__ == "__main__":
    # A positive example SMILES:
    smiles_examples = [
        "C(CCCCCCCC[C@H]([C@H]([C@H](COP(O[C@H]1[C@@H]([C@H]([C@@H]([C@H]([C@H]1O)O)O)O)O)(=O)O)NC(=O)C(CCCCCCCCCCCCCCCCCCCCCCCC)O)O)OCCCCC",
        "CCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(=O)O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@H]1OC1O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]1O)[C@H](O)C(O)CCCCCCCCCCCCCCCC"
    ]
    for smi in smiles_examples:
        result, reason = is_inositol_phosphoceramide(smi)
        print("SMILES:", smi)
        print("Classification:", result)
        print("Reason:", reason)
        print("-----")