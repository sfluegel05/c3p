"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: tripeptide — any oligopeptide consisting of three amino‐acid residues connected by peptide bonds.
This approach works as follows:
  1. The SMILES string is parsed.
  2. We find all occurrences of peptide (amide) bonds using the SMARTS "[C](=O)[N]".
  3. For each match, we unpack the three atoms (carbonyl carbon, oxygen, and amide nitrogen) and discount any bond
     that appears to come from a protecting group (e.g., if the carbonyl carbon has an extra methyl group).
  4. For a linear tripeptide, exactly two “internal” peptide bonds should be present.
  5. In addition, the molecule’s molecular weight and heavy atom count are checked against expected ranges for tripeptides.
If all tests pass, the function returns True along with a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide (three amino acid residues connected by peptide bonds)
    based on its SMILES string.
    
    Heuristic criteria:
      - The molecule must contain exactly two internal peptide (amide) bonds after discounting those which
        come from common protecting groups (such as an acetyl group at the N-terminus).
      - The overall molecular weight should be in a typical tripeptide range (roughly 150-600 Da).
      - The total heavy atom count must also fall within a reasonable range (roughly 15-50).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as a tripeptide,
                     and the second element provides a reason; otherwise False with an explanation.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Calculate molecular weight and heavy atom count as sanity checks.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    # Typical tripeptides have molecular weights between ~150 and 600 Da.
    if not (150 <= mw <= 600):
        return False, f"Molecular weight ({mw:.1f} Da) is outside tripeptide range (150-600 Da)"
    # Typical heavy atom counts for tripeptides fall roughly between 15 and 50.
    if not (15 <= heavy_atoms <= 50):
        return False, f"Heavy atom count ({heavy_atoms}) is outside expected range for tripeptides (15-50)"
        
    # Define a SMARTS for a generic peptide (amide) bond: [C](=O)[N]
    peptide_bond_smarts = "[C](=O)[N]"
    peptide_bond = Chem.MolFromSmarts(peptide_bond_smarts)
    if peptide_bond is None:
        return False, "Could not compile the amide-bond SMARTS pattern"
    
    matches = mol.GetSubstructMatches(peptide_bond)
    internal_bond_count = 0
    
    # For each match, unpack the three atoms: carbonyl carbon, oxygen, and amide nitrogen.
    for match in matches:
        # Unpack carbon index, oxygen index, and nitrogen index.
        c_idx, o_idx, n_idx = match
        c_atom = mol.GetAtomWithIdx(c_idx)
        discount = False
        
        # Check neighbors of the carbonyl carbon to ensure it is not part of a protecting group.
        for nbr in c_atom.GetNeighbors():
            # Skip the neighbor that is the oxygen in the carbonyl group.
            if nbr.GetIdx() == o_idx:
                continue
            # Get the bond between c_atom and the neighbor.
            bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
            if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                continue
            # If a neighbor is a carbon and appears to be a methyl (CH3) group,
            # then it is likely part of an N-terminal protecting (acetyl) group.
            if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() == 3:
                discount = True
                break
        if not discount:
            internal_bond_count += 1
            
    # For a linear tripeptide we expect exactly 2 internal peptide bonds.
    if internal_bond_count != 2:
        return False, f"Found {internal_bond_count} internal peptide bond(s); expected 2 for a tripeptide"

    return True, "Found 2 internal peptide bonds with appropriate molecular size and composition for a tripeptide"

# (Optional) Example usage:
if __name__ == '__main__':
    test_smiles = {
        "trofinetide": "C[C@]1(CCCN1C(=O)CN)C(=O)N[C@@H](CCC(O)=O)C(O)=O",
        "Glu-Glu-Glu": "C(=O)([C@@H](N)CCC(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)CCC(=O)O)CCC(=O)O",
        "N-acetyl-L-tyrosylglycylglycine": "CC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(O)=O",
        "Leu-Thr-Ala": "CC(C)C[C@H](N)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](C)C(O)=O"
    }
    for name, s in test_smiles.items():
        result, reason = is_tripeptide(s)
        print(f"{name}: {result}, {reason}")