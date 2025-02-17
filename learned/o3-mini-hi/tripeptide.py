"""
Classifies: CHEBI:47923 tripeptide
"""
#!/usr/bin/env python
"""
Classifies: tripeptide – any oligopeptide consisting of three amino‐acid residues connected by peptide linkages.
Improved heuristic approach:
  1. Parse the SMILES and check overall molecular weight and heavy atom count.
  2. Identify amide (peptide) bonds using the SMARTS pattern "[C](=O)[N]".
  3. Discount bonds that appear to be part of an N-terminal protection group (heuristic: a CH3 group attached to the carbonyl carbon).
  4. Require exactly 2 valid peptide bonds.
  5. Fragment the molecule on those 2 bonds. In a linear tripeptide the fragmentation should yield exactly 3 fragments.
  
Note: This is a heuristic method and will occasionally mis‐classify molecules.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide (three amino acid residues connected by peptide bonds)
    based on its SMILES string.
    
    Heuristic criteria:
      - Molecular weight and heavy atom count are within generous, expected bounds.
      - Exactly 2 valid peptide (amide) bonds are found using the substructure query "[C](=O)[N]",
        discounting those bonds where the carbonyl carbon is attached to a potential protecting group (e.g. CH3).
      - Fragmentation on these 2 bonds yields exactly 3 fragments, as expected in a linear tripeptide.
    
    Args:
        smiles (str): A SMILES string of the molecule.
        
    Returns:
        (bool, str): A tuple. The first element is True if the molecule is classified as a tripeptide;
                     otherwise False. The second element explains the reasoning.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular size: allow MW in range 150-1000 Da and 15-100 heavy atoms.
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
    
    # Find all substructure matches for the amide bond.
    # Each match is a tuple (c_idx, o_idx, n_idx) where:
    #   c_idx: carbonyl carbon; o_idx: carbonyl oxygen; n_idx: amide nitrogen.
    matches = mol.GetSubstructMatches(peptide_bond_query)
    valid_bonds = []  # Will store tuples (c_idx, n_idx) for bonds we consider valid
    
    for match in matches:
        if len(match) < 3:
            continue  # Expecting three atoms in the match
        c_idx, o_idx, n_idx = match[0], match[1], match[2]
        c_atom = mol.GetAtomWithIdx(c_idx)
        # Discount bonds likely coming from a protecting group:
        # Heuristic: if the carbonyl carbon has a neighbor (other than the carbonyl oxygen)
        # that is a methyl group (a carbon with three attached hydrogens), discount this bond.
        discount = False
        for nbr in c_atom.GetNeighbors():
            # Skip the carbonyl oxygen
            if nbr.GetIdx() == o_idx:
                continue
            if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() == 3:
                discount = True
                break
        if not discount:
            valid_bonds.append((c_idx, n_idx))
    
    if len(valid_bonds) != 2:
        return False, f"Found {len(valid_bonds)} internal peptide bond(s); expected 2 for a tripeptide"
    
    # Next, verify that these two peptide bonds belong to a continuous (linear) chain by performing fragmentation.
    bond_indices = []
    for (c_idx, n_idx) in valid_bonds:
        bond = mol.GetBondBetweenAtoms(c_idx, n_idx)
        if bond is None:
            return False, "Could not retrieve bond between peptide bond atoms"
        bond_indices.append(bond.GetIdx())
    
    # Break the molecule on the 2 peptide bonds. This inserts dummy atoms at the break points.
    frag_mol = Chem.FragmentOnBonds(mol, bond_indices, addDummies=True)
    fragments = Chem.GetMolFrags(frag_mol, asMols=True)
    
    # In a linear tripeptide, breaking the two internal bonds should yield exactly 3 fragments.
    if len(fragments) != 3:
        return False, f"Fragmentation produced {len(fragments)} fragments; peptide bonds may not be connected in a linear chain"
    
    return True, "Valid tripeptide: exactly 2 peptide bonds and fragmentation yields 3 connected fragments"

# (Optional) Example usage and testing.
if __name__ == '__main__':
    test_examples = {
        "trofinetide": "C[C@]1(CCCN1C(=O)CN)C(=O)N[C@@H](CCC(O)=O)C(O)=O",
        "Glu-Glu-Glu": "C(=O)([C@@H](N)CCC(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)CCC(=O)O)CCC(=O)O",
        "N-acetyl-L-tyrosylglycylglycine": "CC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(O)=O",
        "Leu-Thr-Ala": "CC(C)C[C@H](N)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](C)C(O)=O",
        # Additional examples (including some false positives/negatives from evaluation)
        "biotin-valyl-alanyl-aspartyl-fluoromethyl ketone": "[H][C@]12CS[C@@H](CCCCC(=O)N[C@@H](C(C)C)C(=O)N[C@@H](C)C(=O)N[C@@H](CC(O)=O)C(=O)CF)[C@@]1([H])NC(=O)N2",
        "glutathione sulfonate": "S(S(O)(=O)=O)C[C@H](N)C(=O)N(C(=O)CC[C@H](N)C(O)=O)CC(O)=O"
    }
    for name, smi in test_examples.items():
        result, reason = is_tripeptide(smi)
        print(f"{name}: {result} – {reason}")