"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: tripeptide – any oligopeptide consisting of three amino‐acid residues connected by peptide bonds.
Heuristic approach:
  1. The SMILES string is parsed.
  2. Overall molecular size (molecular weight and heavy atom count) is checked (with more generous bounds).
  3. We search for amide (peptide) bonds using the SMARTS “[C](=O)[N]”.
  4. For each match we “discount” those in which the carbonyl carbon appears part of an N‐terminal protecting group 
     (e.g. if a methyl group is attached besides the carbonyl oxygen).
  5. We require that exactly 2 “internal” peptide bonds remain.
  6. Finally we check that these two bonds appear “connected” in that one of the atoms in one match is (through one bond)
     adjacent to an atom that participates in the other match. (In a linear tripeptide, the residue in the middle is “bridged.”
     This is only a rough check.)
If these conditions are met, we classify the structure as a tripeptide.
  
WARNING: This is a heuristic method. Some modified tripeptides (or non‐peptides with two amide bonds)
may be mis‐classified.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide (three amino acid residues connected by peptide bonds)
    based on its SMILES string.
    
    Heuristic criteria:
      - The molecule’s molecular weight is within a relatively generous range (150-1000 Da)
      - The total heavy atom count is within expected bounds (15-100)
      - The molecule contains exactly two (internal) peptide (amide) bonds after discounting those likely 
        arising from protecting groups, as judged by a simple local pattern near the carbonyl carbon.
      - A rudimentary connectivity check is applied: in a linear tripeptide the two peptide bonds should
        “bridge” a middle residue.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as a tripeptide,
                     and the second element provides a reason; otherwise False with an explanation.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular size -- be generous to cover modified residues.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    if not (150 <= mw <= 1000):
        return False, f"Molecular weight ({mw:.1f} Da) is outside acceptable range (150-1000 Da) for modified tripeptides"
    if not (15 <= heavy_atoms <= 100):
        return False, f"Heavy atom count ({heavy_atoms}) is outside expected range (15-100) for tripeptides"
    
    # Define a SMARTS pattern for a generic peptide (amide) bond.
    peptide_bond_smarts = "[C](=O)[N]"
    peptide_bond_query = Chem.MolFromSmarts(peptide_bond_smarts)
    if peptide_bond_query is None:
        return False, "Could not compile the peptide-bond SMARTS pattern"
    
    # Find all substructural matches for the amide bond. Each match is a tuple of atom indices:
    # (carbonyl carbon, carbonyl oxygen, amide nitrogen)
    matches = mol.GetSubstructMatches(peptide_bond_query)
    valid_bonds = []
    
    for match in matches:
        c_idx, o_idx, n_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        
        # Discount bonds that are likely part of an N-terminal blocking/protecting group.
        # Heuristic: if the carbonyl carbon (c_atom) has a non-carbonyl neighbor (other than the oxygen)
        # that is a carbon with three (explicit) hydrogens (i.e. a CH3 group), then discount.
        discount = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetIdx() == o_idx:
                continue  # skip the carbonyl oxygen
            # If neighbor is carbon and appears as a methyl group, discount this bond.
            # (This is only a crude method to remove, e.g., acetyl groups on the peptide N-terminus.)
            if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() == 3:
                discount = True
                break
        if not discount:
            valid_bonds.append((c_idx, n_idx))
    
    if len(valid_bonds) != 2:
        return False, f"Found {len(valid_bonds)} internal peptide bond(s); expected 2 for a tripeptide"
    
    # Optional connectivity check:
    # In a linear tripeptide, the two peptide bonds should belong to one continuous chain.
    # We require that one end (either the amide nitrogen of the first bond or the carbonyl carbon of the second bond)
    # is adjacent (via at least one bond) to an atom that is directly attached to the other bond.
    #
    # Here we simply check that one of the atoms in the first bond (c1 or n1) is within one bond of one of the atoms in the second bond.
    (c1, n1) = valid_bonds[0]
    (c2, n2) = valid_bonds[1]
    
    def is_adjacent(idx1, idx2):
        atom1 = mol.GetAtomWithIdx(idx1)
        for nbr in atom1.GetNeighbors():
            if nbr.GetIdx() == idx2:
                return True
        return False
    
    connectivity_found = (is_adjacent(n1, c2) or is_adjacent(c1, n2) or
                          is_adjacent(n1, n2) or is_adjacent(c1, c2))
    if not connectivity_found:
        return False, "Peptide bonds do not appear to be connected in a continuous chain"
    
    return True, "Found 2 internal peptide bonds with appropriate molecular size and connectivity for a tripeptide"

# (Optional) Example usage and testing:
if __name__ == '__main__':
    test_smiles = {
        "trofinetide": "C[C@]1(CCCN1C(=O)CN)C(=O)N[C@@H](CCC(O)=O)C(O)=O",
        "Glu-Glu-Glu": "C(=O)([C@@H](N)CCC(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)CCC(=O)O)CCC(=O)O",
        "N-acetyl-L-tyrosylglycylglycine": "CC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(O)=O",
        "Leu-Thr-Ala": "CC(C)C[C@H](N)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](C)C(O)=O",
        # A few other examples (including some known false positives/negatives from evaluation)
        "LSM-7408": "C[C@H]1CN(C(=O)C2=CC(=CN=C2O[C@@H]1CN(C)C(=O)CC3=CN=CC=C3)C=CC4=CC=CC=C4)[C@@H](C)CO",
        "biotin-valyl-alanyl-aspartyl-fluoromethyl ketone": "[H][C@]12CS[C@@H](CCCCC(=O)N[C@@H](C(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](CC(O)=O)C(=O)CF)[C@@]1([H])NC(=O)N2"
    }
    for name, s in test_smiles.items():
        result, reason = is_tripeptide(s)
        print(f"{name}: {result} – {reason}")