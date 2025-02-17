"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: tripeptide — any oligopeptide consisting of three amino‐acid residues connected by peptide linkages.

This improved approach works as follows:
  1. The SMILES string is parsed.
  2. We find all occurrences of amide bonds using the SMARTS "[C](=O)[N]".
  3. For each amide bond found we discount those that appear to come from an N‐terminal protecting group 
     (e.g. an acetyl group is often –CH3C(=O)–, so if the carbonyl carbon is bonded to a methyl group aside 
     from the amide nitrogen then we discount that match).
  4. For a linear tripeptide exactly two “internal” peptide bonds (after discounting) should be present.
  5. In addition the molecule’s molecular weight and heavy atom count are checked against ranges expected for tripeptides.
  
If all tests pass, the function returns True along with a reason. Otherwise it returns False with an explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide (three amino acid residues connected by peptide bonds)
    based on its SMILES string.
    
    Heuristic criteria:
      - The molecule must contain exactly two internal peptide (amide) bonds.
      - We use a SMARTS pattern "[C](=O)[N]" to find amide bonds while discounting those coming from protecting groups,
        for example an acetyl group at the N‐terminus.
      - The overall molecular weight must be in a tripeptide range (roughly 150-600 Da).
      - The total number of non‐hydrogen atoms (heavy atoms) must be within a reasonable range (roughly 15-50).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as a tripeptide,
                     and the second element provides a reason for the decision.
                     If the SMILES is invalid or criteria are not met, returns False with an explanation.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Get molecular weight and heavy atom count; these will help filter out false matches.
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    # Expect typical tripeptides to have molecular weights between ~150 and 600 Da
    if not (150 <= mw <= 600):
        return False, f"Molecular weight ({mw:.1f} Da) is outside tripeptide range (150-600 Da)"
    # Expect typical tripeptides to have roughly 15 to 50 heavy atoms (adjustable)
    if not (15 <= heavy_atoms <= 50):
        return False, f"Heavy atom count ({heavy_atoms}) is outside expected range for tripeptides (15-50)"
        
    # Define a SMARTS for a generic amide bond: [C](=O)[N]
    peptide_bond_smarts = "[C](=O)[N]"
    peptide_bond = Chem.MolFromSmarts(peptide_bond_smarts)
    if peptide_bond is None:
        return False, "Could not compile the amide-bond SMARTS pattern"
    
    matches = mol.GetSubstructMatches(peptide_bond)
    # For each match, check if the carbonyl carbon seems to be part of a protecting group.
    # In particular, if the carbonyl carbon (first atom in the match) has a neighbor (other than the amide N)
    # that is a methyl group (CH3) then discount that match.
    internal_bond_count = 0
    for match in matches:
        c_idx, n_idx = match  # carbonyl carbon and its bonded nitrogen
        c_atom = mol.GetAtomWithIdx(c_idx)
        discount = False
        # Loop over neighbors of the carbonyl carbon
        for nbr in c_atom.GetNeighbors():
            # Skip the neighbor that is the amide nitrogen in this bond
            if nbr.GetIdx() == n_idx:
                continue
            # Also skip neighbors connected by a double bond to oxygen (the carbonyl O)
            bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), nbr.GetIdx())
            if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                continue
            # If a neighbor is carbon and looks like a CH3 group (has 3 total H's), discount this bond;
            # this is to handle the common N-terminal acetyl group.
            if nbr.GetAtomicNum() == 6:
                # Get the total number of attached hydrogens
                if nbr.GetTotalNumHs() == 3:
                    discount = True
                    break
        if not discount:
            internal_bond_count += 1
            
    if internal_bond_count != 2:
        return False, f"Found {internal_bond_count} internal peptide bond(s); expected 2 for a tripeptide"

    # If all tests pass, we have a tripeptide.
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