"""
Classifies: CHEBI:47923 tripeptide
"""
"""
Classifies: tripeptide — any oligopeptide consisting of three amino‐acid residues connected by peptide linkages.
The function is_tripeptide takes a SMILES string and returns a boolean value plus a reason for the classification.
"""

from rdkit import Chem

def is_tripeptide(smiles: str):
    """
    Determines if a molecule is a tripeptide (three amino acid residues connected by peptide bonds)
    based on its SMILES string.
    
    Our heuristic approach is as follows:
      1. Parse the SMILES string.
      2. Use a SMARTS for a peptide-bond fragment "[C](=O)[N][C]". 
         (This pattern is intended to pick up an amide bond where the nitrogen is bound to a carbon,
          as expected in the connective backbone.)
      3. Count the number of “internal peptide bonds”. In a linear tripeptide the peptide backbone 
         should have exactly two internal peptide bonds linking three residues.
      4. We attempt to discount peptide bonds coming from N‐terminal modifications (commonly acetyl groups).
         For each match we check the carbonyl carbon (first atom in the match) – if it is attached 
         to a CH3 group (which is typical for an acetyl substituent), we do not count that peptide bond.
      5. If the count of internal peptide bonds equals 2, we return True with an explanation;
         otherwise we return False with a reason.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as a tripeptide,
                     and the second element provides a brief explanation.
                     If the SMILES is invalid or the molecule fails the tripeptide heuristic,
                     returns False along with a reason.
    """
    # Parse the SMILES string using rdkit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Define a SMARTS pattern for a generic amide bond found in peptide backbones.
    # This pattern matches a carbon (the carbonyl carbon) that is double-bonded to oxygen,
    # then single-bonded to a nitrogen which in turn is bonded to a carbon.
    peptide_bond_smarts = "[C](=O)[N][C]"
    peptide_bond = Chem.MolFromSmarts(peptide_bond_smarts)
    if peptide_bond is None:
        return False, "Could not compile peptide bond SMARTS pattern"
        
    # Find all matches (each match is a tuple of atom indices corresponding to [C, O?, N, C]).
    matches = mol.GetSubstructMatches(peptide_bond)
    
    # We'll count only those peptide bonds that are likely to be “internal” in the backbone.
    internal_count = 0
    for match in matches:
        # In our match the atoms are:
        # match[0]: the carbonyl carbon; match[1]: the amide nitrogen; match[2]: the following carbon.
        # (Note: our SMARTS implies a bond C(=O)-N-C)
        c_carbon = mol.GetAtomWithIdx(match[0])
        # We try to discount amide bonds that originate from a protecting group.
        # For example, an acetyl group is CH3-C(=O)-N.... 
        # We check the neighbors of the carbonyl carbon (excluding the oxygen from the carbonyl).
        discount = False
        for nbr in c_carbon.GetNeighbors():
            # Exclude neighbors connected by a double bond to oxygen (the carbonyl O)
            bond = mol.GetBondBetweenAtoms(c_carbon.GetIdx(), nbr.GetIdx())
            if bond is not None and bond.GetBondTypeAsDouble() == 2.0:
                continue
            # If a neighbor has atomic number 6 (carbon) and 
            # (typically in an acetyl group it is a CH3 — often with low degree),
            # we discount this match as coming from an N-terminal protecting group.
            if nbr.GetAtomicNum() == 6 and nbr.GetDegree() <= 2:
                discount = True
                break
        if not discount:
            internal_count += 1

    # For a linear tripeptide there should be exactly two internal peptide bonds connecting the three residues.
    if internal_count == 2:
        return True, "Found 2 backbone peptide bonds connecting three amino acid residues"
    else:
        return False, f"Found {internal_count} internal peptide bond(s); expected 2 for a tripeptide"

# (Optional) Example usage:
if __name__ == '__main__':
    # List a few SMILES strings for testing (some from the examples provided)
    test_smiles = {
        "trofinetide": "C[C@]1(CCCN1C(=O)CN)C(=O)N[C@@H](CCC(O)=O)C(O)=O",
        "Glu-Glu-Glu": "C(=O)([C@@H](N)CCC(=O)O)N[C@H](C(=O)N[C@H](C(=O)O)CCC(=O)O)CCC(=O)O",
        "N-acetyl-L-tyrosylglycylglycine": "CC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)NCC(=O)NCC(O)=O",
        "Leu-Thr-Ala": "CC(C)C[C@H](N)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](C)C(O)=O"
    }
    for name, s in test_smiles.items():
        result, reason = is_tripeptide(s)
        print(f"{name}: {result}, {reason}")