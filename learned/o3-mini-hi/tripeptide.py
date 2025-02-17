"""
Classifies: CHEBI:47923 tripeptide
"""
#!/usr/bin/env python
"""
Classifies: tripeptide – any oligopeptide consisting of three amino‐acid residues connected by peptide linkages.
An amino acid residue is assumed to be connected through “main‐chain” peptide bonds.
This improved heuristic works as follows:
  1. Parse the molecule and check overall molecular weight and heavy atom count.
  2. Find all candidate peptide bonds via the SMARTS pattern "[C](=O)[N]".
  3. (A simple discount is applied to ignore bonds that seem to come from a protecting group:
     if the carbonyl carbon is directly attached to a methyl group (CH3) other than the amide oxygen, we discount it.)
  4. Since extra amide bonds may be present, try all combinations of 2 candidate bonds.
  5. For each combination, fragment the molecule on those bonds. In a linear tripeptide the fragmentation yields exactly 3 fragments.
     Moreover, dummy atoms (atomic num=0) are inserted at the break points so that one fragment should have 2 dummy atoms (central residue)
     and the other two fragments exactly one dummy atom each.
  6. If any pair of bonds meets these criteria we return True, otherwise False.
  
Note: This heuristic method is not perfect and is sensitive to protecting groups and non‐canonical linkages.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import itertools

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide (three amino acid residues connected by peptide bonds)
    based on its SMILES string.
    
    Heuristic criteria:
      - The molecular weight is between 200 and 1000 Da and the heavy atom count is between 15 and 100.
      - Candidate amide (peptide) bonds are found using the SMARTS pattern "[C](=O)[N]". Bonds are discounted
        if the carbonyl carbon has a nonamide neighbor that is a methyl group (CH3), which may indicate an N‑terminal
        protection group.
      - Of all candidate bonds, we try all combinations of two bonds. If fragmenting on those bonds yields exactly 3 fragments,
        and the dummy-atom counts in the fragments are exactly [1, 1, 2] (i.e. the middle residue has two break points)
        then we call that a valid, linear tripeptide.
    
    Args:
        smiles (str): A SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple. The first element is True if the molecule is classified as a tripeptide;
                     otherwise False. The second element explains the reasoning.
    """
    # First, parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Overall molecular size check. Most tripeptides fall roughly into this range.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    if not (200 <= mw <= 1000):
        return False, f"Molecular weight ({mw:.1f} Da) is outside acceptable range (200-1000 Da) for tripeptides"
    if not (15 <= heavy_atoms <= 100):
        return False, f"Heavy atom count ({heavy_atoms}) is outside expected range (15-100) for tripeptides"
    
    # Define a SMARTS pattern for a generic peptide (amide) bond.
    peptide_bond_smarts = "[C](=O)[N]"
    peptide_bond_query = Chem.MolFromSmarts(peptide_bond_smarts)
    if peptide_bond_query is None:
        return False, "Could not compile the peptide-bond SMARTS pattern"
    
    # Find all substructure matches for the amide bond.
    # Each match should include at least three atoms: carbonyl carbon, carbonyl oxygen, and amide nitrogen.
    all_matches = mol.GetSubstructMatches(peptide_bond_query)
    candidate_bonds = []
    for match in all_matches:
        if len(match) < 3:
            continue
        c_idx, o_idx, n_idx = match[0], match[1], match[2]
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Discount candidate if the carbonyl carbon has a neighbor (other than the C=O oxygen)
        # with atomic number 6 (carbon) and three attached hydrogens (a CH3 group),
        # which often indicates a protecting (acetyl) group.
        discount = False
        for nbr in c_atom.GetNeighbors():
            if nbr.GetIdx() == o_idx:
                continue
            if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() == 3:
                discount = True
                break
        if not discount:
            candidate_bonds.append( (c_idx, n_idx) )
    
    if len(candidate_bonds) < 2:
        return False, f"Found only {len(candidate_bonds)} candidate peptide bond(s); need at least 2 for a tripeptide"
    
    # Now try every combination of 2 candidate bonds.
    # We will fragment on the corresponding bonds and check if we get exactly 3 fragments
    # with the expected pattern of dummy atoms: two fragments having one dummy each (terminals) and one with two dummies (center).
    valid_linear_chain_found = False
    reason_attempts = []
    for bond_pair in itertools.combinations(candidate_bonds, 2):
        bond_indices = []
        skip_pair = False
        for (c_idx, n_idx) in bond_pair:
            bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
            if bond is None:
                skip_pair = True
                break
            bond_indices.append(bond.GetIdx())
        if skip_pair:
            continue
        
        # Perform fragmentation on these two bonds.
        frag_mol = Chem.FragmentOnBonds(mol, bond_indices, addDummies=True)
        fragments = Chem.GetMolFrags(frag_mol, asMols=True)
        if len(fragments) != 3:
            reason_attempts.append(f"Fragmentation on bonds {bond_indices} yielded {len(fragments)} fragments")
            continue
        
        # Count dummy atoms (dummy atoms have atomic number 0) in each fragment.
        dummy_counts = [sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 0) for frag in fragments]
        # In a linear backbone: the two terminal residues get one dummy atom each; the central residue gets two.
        if sorted(dummy_counts) == [1, 1, 2]:
            valid_linear_chain_found = True
            break
        else:
            reason_attempts.append(f"Dummy atoms per fragment were {dummy_counts} for bonds {bond_indices}")
    
    if valid_linear_chain_found:
        return True, "Valid tripeptide: a pair of peptide bonds yields 3 fragments with dummy distribution [1,1,2]"
    else:
        detailed = "; ".join(reason_attempts[-3:])  # include up to last few reasons
        return False, f"No candidate peptide bond pair produced a valid linear tripeptide backbone. {detailed}"

# (Optional) Example usage and testing.
if __name__ == '__main__':
    # Some examples from the provided list:
    examples = {
        "trofinetide": "C[C@]1(CCCN1C(=O)CN)C(=O)N[C@@H](CCC(O)=O)C(O)=O",
        "Glu-Glu-Glu": "C(=O)([C@@H](N)CCC(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)CCC(=O)O)CCC(=O)O",
        "Leu-Thr-Ala": "CC(C)C[C@H](N)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](C)C(O)=O",
        "N-acetyl-L-tyrosylglycylglycine": "CC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(O)=O",
        "His-Gly-Ser (false positive case?)": "O=C(NCC(=O)N[C@@H](CO)C(O)=O)[C@@H](N)CC=1NC=NC1",
        "Ac-D-Glu-His-D-Pro-NH2 (false negative case?)": "CC(=O)N[C@H](CCC(O)=O)C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N1CCC[C@@H]1C(N)=O"
    }
    
    for name, smi in examples.items():
        res, reason = is_tripeptide(smi)
        print(f"{name}: {res} – {reason}")