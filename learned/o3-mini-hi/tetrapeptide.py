"""
Classifies: CHEBI:48030 tetrapeptide
"""
#!/usr/bin/env python3
"""
Classifies: Tetrapeptide
Definition: Any molecule that contains four amino-acid residues connected by peptide linkages.
Heuristic:
  - A linear tetrapeptide should have three peptide bonds linking four residues.
  - Many tetrapeptides will contain four alpha-carbons bearing an amino group and adjacent to a carbonyl.
This implementation first searches for backbone peptide bond fragments using the SMARTS
"C(=O)N[C@H]" and "C(=O)N[C@@H]". It then also searches for “alpha-residue” fragments defined by
[C@H](N)C(=O) or [C@@H](N)C(=O). If either the peptide-bond count is exactly 3 or there are exactly 4
alpha-residues, and additional properties (molecular weight and nitrogen count) fall into a typical range,
the molecule is classified as a tetrapeptide.
Note:
  This heuristic is imperfect and may mis‐classify modified or cyclic peptides as well as non-peptidic molecules
  that accidentally contain a similar motif.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tetrapeptide(smiles: str):
    """
    Determines if a molecule is a tetrapeptide based on its SMILES string.
    The algorithm uses two complementary criteria:
      1. Count peptide bonds by matching SMARTS for backbone fragments.
         A linear tetrapeptide is expected to have three peptide bonds linking four alpha-carbons.
      2. Count “alpha‐residues” using SMARTS for the core of an amino acid residue, i.e.
         a chiral carbon attached to an amino group and a carbonyl (which should be four in a tetrapeptide).
    Additionally, properties like molecular weight and the number of nitrogen atoms are checked.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is likely a tetrapeptide, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --------------------------
    # Criterion 1: Count peptide bonds in the backbone.
    # We use two SMARTS patterns to capture the fragment "C(=O)N[C@H]" and "C(=O)N[C@@H]".
    pat1 = Chem.MolFromSmarts("C(=O)N[C@H]")
    pat2 = Chem.MolFromSmarts("C(=O)N[C@@H]")
    if pat1 is None or pat2 is None:
        return False, "Error creating backbone SMARTS patterns"
    
    peptide_bonds = set()
    # RDKit SMARTS for "C(=O)N[C@H]" will typically match four atoms:
    # carbonyl C, carbonyl O (implied in the =O), amide N and the following chiral carbon.
    for match in mol.GetSubstructMatches(pat1):
        # Ensure match length is as expected.
        if len(match) >= 3:
            c_idx = match[0]  # carbonyl carbon
            n_idx = match[2]  # amide nitrogen
            peptide_bonds.add((c_idx, n_idx))
    for match in mol.GetSubstructMatches(pat2):
        if len(match) >= 3:
            c_idx = match[0]
            n_idx = match[2]
            peptide_bonds.add((c_idx, n_idx))
    n_peptide = len(peptide_bonds)
    
    # --------------------------
    # Criterion 2: Count alpha-residue cores.
    # "Alpha-residue" here is defined as a chiral carbon attached to an amino (N) and a carbonyl C.
    # We use SMARTS "[C@H](N)C(=O)" and "[C@@H](N)C(=O)".
    alpha_pat1 = Chem.MolFromSmarts("[C@H](N)C(=O)")
    alpha_pat2 = Chem.MolFromSmarts("[C@@H](N)C(=O)")
    if alpha_pat1 is None or alpha_pat2 is None:
        return False, "Error creating alpha-residue SMARTS patterns"
    
    alpha_residues = set()
    for match in mol.GetSubstructMatches(alpha_pat1):
        # match returns indices for the chiral carbon, the attached N, and the carbonyl carbon.
        if len(match) == 3:
            alpha_residues.add(match[0])
    for match in mol.GetSubstructMatches(alpha_pat2):
        if len(match) == 3:
            alpha_residues.add(match[0])
    n_alpha = len(alpha_residues)
    
    # --------------------------
    # Additional property checks:
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    # Typical tetrapeptides are in the range 250-1500 Da. (Modified peptides might lie outside this range.)
    if not (250 <= mol_wt <= 1500):
        return False, "Molecular weight {:.1f} Da out of typical tetrapeptide range (250-1500 Da).".format(mol_wt)
    
    # Count nitrogen atoms, which should be within a plausible range
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if not (4 <= n_nitrogen <= 20):
        return False, "Found {} nitrogen atom(s), which is not typical for a tetrapeptide.".format(n_nitrogen)
    
    # --------------------------
    # Combine the criteria:
    # For a linear tetrapeptide we expect:
    #   - exactly 3 backbone peptide bonds (criterion 1)
    #   OR
    #   - exactly 4 distinct alpha-residues (criterion 2)
    if n_peptide == 3:
        return True, ("Found {} peptide bond(s) matching the backbone motif and overall properties "
                      "consistent with a tetrapeptide.".format(n_peptide))
    elif n_alpha == 4:
        return True, ("Found {} alpha-residue fragment(s) consistent with 4 amino-acid residues "
                      "and overall properties consistent with a tetrapeptide.".format(n_alpha))
    else:
        return False, ("Found {} peptide bond(s) and {} alpha-residue fragment(s); expected either 3 peptide bonds "
                       "or 4 alpha-residues for a linear tetrapeptide.".format(n_peptide, n_alpha))

# For quick testing, uncomment one or more test cases below:
if __name__ == '__main__':
    # Example test: Lys-Thr-Trp-Tyr (a true tetrapeptide)
    test_smiles = "C[C@@H](O)[C@H](NC(=O)[C@@H](N)CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccc(O)cc1)C(O)=O"
    result, reason = is_tetrapeptide(test_smiles)
    print(result, reason)