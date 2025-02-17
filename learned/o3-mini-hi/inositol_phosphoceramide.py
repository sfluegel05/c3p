"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: Inositol Phosphoceramide
Definition: A phosphosphingolipid in which an inositol residue and the ceramide moiety are linked via a phosphodiester bridge.
The ceramide moiety contains substituents (typically an amide bond and long aliphatic chains).
In this implementation:
  - A strict SMARTS pattern is used to capture a cyclohexane bearing an OH on every carbon.
  - We locate the unique inositol ring and record its atom indices.
  - We collect all matches for an amide bond (C(=O)N) and record the carbonyl carbons.
  - For every phosphorus atom, we inspect its bonds: it must have at least one double‐bonded oxygen (P=O) and two single-bonded oxygen neighbors.
    One of these oxygens must connect to an inositol atom and another must be bonded to one of the amide carbonyl carbons.
  - We also require the overall molecule to have a minimum molecular weight and number of carbons.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    The method checks:
      1) The presence of an inositol ring substructure: a cyclohexane that carries an OH on each carbon.
      2) A phosphodiester linkage in which a phosphorus atom is connected via oxygen to both the inositol and
         to a ceramide moiety (the latter being indicated by the presence of an amide bond, C(=O)N).
      3) Sufficient overall molecular size.
    
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
    # Use a SMARTS for myo-inositol: a cyclohexane with an OH on every carbon.
    inositol_smarts = "OC1C(O)C(O)C(O)C(O)C1O"
    inositol_query = Chem.MolFromSmarts(inositol_smarts)
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        return False, "Missing inositol ring substructure"
        
    # For further connectivity check, gather all atom indices that belong to any inositol match.
    inositol_atom_indices = set()
    for match in inositol_matches:
        inositol_atom_indices.update(match)
    
    # --- Step 2: Find amide bonds (C(=O)N) indicative of the ceramide part ---
    amide_smarts = "C(=O)N"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide_query):
        return False, "Missing ceramide amide bond"
    
    # Get set of carbon atoms that appear as the carbonyl carbon in any amide match.
    amide_matches = mol.GetSubstructMatches(amide_query)
    carbonyl_carbons = set(match[0] for match in amide_matches)  # first atom is C
    
    # --- Step 3: Check phosphate linkage ---
    # We iterate over phosphorus atoms and check if one of its single-bonded O neighbors leads to 
    # an inositol ring and another to a carbonyl carbon from an amide bond.
    phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue  # not phosphorus
        p_atom = atom
        # Get bonds from phosphorus
        bonds = p_atom.GetBonds()
        double_bonded_os = []
        single_bonded_os = []
        for bond in bonds:
            # Look only at bonds to oxygen atoms
            nbr = bond.GetOtherAtom(p_atom)
            if nbr.GetAtomicNum() != 8:
                continue
            bt = bond.GetBondType()
            if bt == Chem.rdchem.BondType.DOUBLE:
                double_bonded_os.append(nbr)
            elif bt == Chem.rdchem.BondType.SINGLE:
                single_bonded_os.append(nbr)
        
        # We require at least one double-bonded oxygen (P=O) and at least two single bonds.
        if len(double_bonded_os) < 1 or len(single_bonded_os) < 2:
            continue
        
        inositol_linked = False
        ceramide_linked = False
        # For each single-bonded oxygen, inspect its other neighbor(s)
        for oxy in single_bonded_os:
            for nbr in oxy.GetNeighbors():
                if nbr.GetIdx() == p_atom.GetIdx():
                    continue
                # If the neighbor is part of the inositol ring, mark
                if nbr.GetIdx() in inositol_atom_indices:
                    inositol_linked = True
                # If the neighbor is a carbon and is part of an amide bond, mark
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in carbonyl_carbons:
                    ceramide_linked = True
            # When both connections are found, this phosphate qualifies.
            if inositol_linked and ceramide_linked:
                phosphate_found = True
                break
        if phosphate_found:
            break
    
    if not phosphate_found:
        return False, "Phosphodiester linkage (connecting phosphate to both inositol and ceramide) not detected"
    
    # --- Step 4: Check overall molecular size ---
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for inositol phosphoceramide"
    
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Too few carbons ({carbon_count}) to be an inositol phosphoceramide"
    
    return True, "Molecule contains an inositol ring, the phosphate bridges inositol to a ceramide moiety (via an amide bond), and has sufficient molecular size"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)CCCCCCCCCCCCCC"
    result, reason = is_inositol_phosphoceramide(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)