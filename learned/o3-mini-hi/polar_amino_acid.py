"""
Classifies: CHEBI:26167 polar amino acid
"""
#!/usr/bin/env python
"""
Classifies: Polar Amino Acid
Definition: Any amino acid whose side chain is capable of forming one or more hydrogen bonds.
A polar amino acid is defined as a free (non‐peptide) amino acid whose side chain
contains a polar atom (O or N) or a –SH group from sulfur with attached hydrogen.
The algorithm:
  1. Searches for the α–carbon: a carbon atom having exactly three heavy atom (non‐hydrogen) neighbors.
     These three must include one amino nitrogen (N), one carboxyl carbon (bonded via a double bond to oxygen)
     and one branch that is taken as the side chain. Glycine is excluded (it has no side chain).
  2. If more than one such backbone candidate is found, the molecule is assumed to be peptide-like.
  3. Additionally, we check that the molecular weight is below a threshold (here, 250 Da) as free amino acids 
     are small molecules. Larger molecules are likely peptides.
  4. Finally, we traverse the side chain (DFS from the designated side-chain branch) and also require that
     the side chain fragment has only a limited number of heavy atoms (here, no more than 12) to ensure it is not an extended chain.
  5. If a polar atom (N, O, or a –SH group from S with an H) is found in the side-chain, the molecule is called polar.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polar_amino_acid(smiles: str):
    """
    Determines whether the input molecule (given by a SMILES string) is a polar amino acid.
    
    It uses these criteria:
      1. There is exactly one putative free amino acid backbone (non-glycine) display an α–carbon (C)
         with exactly three heavy-atom neighbors.
      2. Among these neighbors, one must be an amino nitrogen (N), one must be a carboxyl carbon 
         (with at least one double bonded oxygen) and one is the start of the side chain.
      3. The overall molecular weight is small (<250 Da) to prevent peptide fragments from being misclassified.
      4. The side chain (found by DFS from the side-chain branch, avoiding the α–carbon) must not be overly large (<=12 heavy atoms).
      5. The side chain must include at least one polar feature (oxygen, nitrogen, or –SH from sulfur with H).
    
    Args:
         smiles (str): SMILES string representation of the molecule.
         
    Returns:
         (bool, str): Tuple of a boolean (True if it is a polar amino acid, False otherwise)
                      and a string with the reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular weight (free amino acids are small)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt > 250:
        return False, f"Molecular weight ({mol_wt:.1f} Da) too high for a free amino acid"
        
    candidates = []
    
    # Look for candidate α–carbons (a carbon with exactly three heavy neighbors).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # α–carbon must be carbon
            continue
        # Identify heavy-atom (non-hydrogen) neighbors
        neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(neighbors) != 3:  # Should have three heavy neighbors (glycine has 2, so it's excluded)
            continue

        amino_found = None
        carboxyl_found = None
        sidechain_found = None

        # Inspect each neighbor.
        for nbr in neighbors:
            anum = nbr.GetAtomicNum()
            if anum == 7:
                amino_found = nbr
            elif anum == 6:
                # Possibility for carboxyl: check for a double bonded oxygen in its neighbors (excluding back to α–carbon)
                oxy_double = False
                for nn in nbr.GetNeighbors():
                    if nn.GetIdx() == atom.GetIdx():
                        continue
                    if nn.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nbr.GetIdx(), nn.GetIdx())
                        if bond and bond.GetBondType() == Chem.BondType.DOUBLE:
                            oxy_double = True
                            break
                if oxy_double:
                    carboxyl_found = nbr
                else:
                    sidechain_found = nbr
            else:
                # any other heavy atom is assigned to sidechain.
                sidechain_found = nbr

        # Require all three kinds of neighbor for a valid free amino acid backbone.
        if amino_found is None or carboxyl_found is None or sidechain_found is None:
            continue

        candidates.append({
            'alpha': atom,
            'amino': amino_found,
            'carboxyl': carboxyl_found,
            'sidechain': sidechain_found
        })
    
    if len(candidates) == 0:
        return False, "No valid free amino acid backbone (non-glycine) found"
    if len(candidates) > 1:
        return False, "Multiple amino acid backbones detected; molecule appears to be a peptide or larger compound"
    
    # Use the sole candidate found.
    backbone = candidates[0]
    alpha = backbone['alpha']
    side_start = backbone['sidechain']
    
    # Traverse the side chain (avoid going back to the α–carbon).
    side_atoms = set()
    stack = [side_start.GetIdx()]
    while stack:
        curr_idx = stack.pop()
        if curr_idx in side_atoms:
            continue
        side_atoms.add(curr_idx)
        curr_atom = mol.GetAtomWithIdx(curr_idx)
        for nbr in curr_atom.GetNeighbors():
            # Do not go back to the α–carbon.
            if nbr.GetIdx() == alpha.GetIdx():
                continue
            if nbr.GetIdx() not in side_atoms:
                stack.append(nbr.GetIdx())
    
    # For a free amino acid, the side-chain portion should be small.
    # If the number of heavy atoms in the side chain is too high, suspect a peptide derivative.
    side_heavy_count = sum(1 for idx in side_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1)
    if side_heavy_count > 12:
        return False, f"Side chain too large ({side_heavy_count} heavy atoms) for a typical free amino acid"
    
    # Now, check for polar features in the side chain.
    polar_features = []
    for idx in side_atoms:
        atom = mol.GetAtomWithIdx(idx)
        anum = atom.GetAtomicNum()
        # Oxygen and nitrogen are polar.
        if anum in (7, 8):
            polar_features.append(atom.GetSymbol())
        elif anum == 16:  # sulfur: count only if attached to at least one hydrogen (–SH group).
            if atom.GetTotalNumHs() > 0:
                polar_features.append("SH")
                
    if polar_features:
        return True, "Side chain contains polar feature(s): " + ", ".join(polar_features)
    else:
        return False, "Side chain does not contain polar functional group(s) capable of hydrogen bonding"

# If run as a script, example test cases.
if __name__ == "__main__":
    test_cases = [
        ("NC(CCC(N)=O)C(O)=O", "glutamine"),
        ("NC(CO)C(O)=O", "serine"),
        ("N[C@H](CS)C(O)=O", "D-cysteine"),
        ("N[C@@H](Cc1c[nH]cn1)C(O)=O", "L-histidine"),
        ("N[C@@H](CCC(N)=O)C(O)=O", "L-glutamine"),
        ("OC(=O)[C@](N([2H])[2H])(C(C(CN=C(N)N)([2H])[2H])([2H])[2H])[2H]", "L-arginine-d7"),
        ("NCCCC[C@@H](N)C(O)=O", "D-lysine"),
        ("N[C@H](CC(O)=O)C(O)=O", "D-aspartic acid"),
        ("NC(Cc1ccc(O)cc1)C(O)=O", "tyrosine"),
        ("N[C@H](CO)C(O)=O", "D-serine"),
        ("NC(Cc1c[nH]c2ccccc12)C(O)=O", "tryptophan"),
        ("N[C@@H](CCCNC(N)=N)C(O)=O", "L-arginine")
    ]
    
    for smi, name in test_cases:
        result, reason = is_polar_amino_acid(smi)
        print(f"SMILES: {smi}\nNAME: {name}\nResult: {result} ({reason})\n")