"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
"""
Classifies: Short‐chain fatty acid according to the definition:
"An aliphatic monocarboxylic acid with a chain length of less than C6.
If any non‐hydrocarbon substituent is present, the compound is not normally regarded as a short‐chain fatty acid."
In our implementation we require:
  – exactly one carboxylic acid group (SMARTS: "[CX3](=O)[OX2H]")
  – acyclic (no rings)
  – the total number of carbon atoms in the molecule (which equals the acid carbon plus the alkyl chain) <= 6
  – no oxygen atoms appear outside of the carboxyl group (i.e. no additional non‐hydrocarbon substituents)
"""

from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is defined as an aliphatic monocarboxylic acid with a chain length
    (total carbon count including the carboxyl carbon) of less than 7 (i.e. up to 6 carbons)
    and with no non-hydrocarbon substituents outside the carboxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as a short-chain fatty acid, False otherwise.
        str: A message describing the reason for the classification.
    """
    # Parse the SMILES string into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # The molecule should not contain any rings (it must be aliphatic).
    if mol.GetRingInfo().NumRings() != 0:
        return False, "Molecule contains rings; expected an acyclic (aliphatic) acid."
    
    # Define a SMARTS for the carboxylic acid group.
    # This pattern matches a carboxyl group: a carbon bonded to a double-bonded oxygen and a hydroxyl oxygen.
    ca_group_smarts = "[CX3](=O)[OX2H]"
    ca_group = Chem.MolFromSmarts(ca_group_smarts)
    
    # Find matches of the carboxyl group in the molecule.
    ca_matches = mol.GetSubstructMatches(ca_group)
    if not ca_matches:
        return False, "No carboxylic acid functional group found"
    if len(ca_matches) != 1:
        return False, f"Expected exactly one carboxyl group, found {len(ca_matches)}"
    
    # Mark the atoms that are part of the carboxyl group (both the carbon and the two oxygens).
    ca_atoms = set(ca_matches[0])
    
    # Count the total number of carbon atoms in the molecule.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # According to the definition, the chain length (i.e. total carbons including the acid carbon) should be < 7.
    if carbon_count >= 7:
        return False, f"Total carbon count is {carbon_count}, which exceeds allowed short-chain length (max 6)"
    
    # Check that no additional non-hydrocarbon substituents (i.e. heteroatoms) are present.
    # Here we allow only carbon and hydrogen outside of the carboxyl group.
    for atom in mol.GetAtoms():
        # If the atom is in the carboxyl group, skip it.
        if atom.GetIdx() in ca_atoms:
            continue
        # Allow only carbon (atomic number 6) or (implicit) hydrogen.
        # Note: Hydrogens are usually implicit, so we check explicit atoms.
        if atom.GetAtomicNum() != 6:
            return False, f"Found a non-hydrocarbon substituent: atom {atom.GetSymbol()} not allowed"
    
    return True, "Molecule contains one carboxyl group, is acyclic, has ≤6 carbons, and no extraneous substituents."

# Example usage (for testing purposes):
if __name__ == '__main__':
    test_smiles = [
        "CC(C)(C)C(O)=O",  # pivalic acid
        "CCC(O)=O",       # propionic acid
        "CCCC(O)=O",      # butyric acid
        "CCCCC(O)=O",     # valeric acid (5 carbons total, qualifies)
        "CCCC(C)C(O)=O",  # 2-methylvaleric acid (6 carbons total, qualifies)
        "OCCC(O)=O",      # 3-hydroxypropionic acid (has an extra OH -> disqualified)
        "[H][C@@]12[C@H](...)",  # heliosupine-like, will fail due to rings & extra atoms.
    ]
    for s in test_smiles:
        valid, reason = is_short_chain_fatty_acid(s)
        print(f"SMILES: {s} -> {valid}; Reason: {reason}")