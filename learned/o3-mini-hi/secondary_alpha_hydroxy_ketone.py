"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: Secondary Alpha-Hydroxy Ketone (Acyloin)

Definition:
  A secondary α‐hydroxy ketone (acyloin) is defined as a molecule containing the motif:
    R–CH(OH)–C(=O)–R′,
  where the “CH(OH)” center is secondary (has exactly one hydrogen and two heavy-atom neighbors)
  and the carbonyl carbon is doubly bonded to an oxygen and is not part of a free acid.

The strategy is:
  1. Use a SMARTS pattern with atom mapping to match the motif.
  2. From each match (which now returns 5 atoms), extract the alcohol (mapped as :alc) and ketone carbon (mapped as :co).
  3. Verify that the alcohol carbon carries exactly one hydrogen.
  4. Check that the ketone carbon has a proper double-bonded oxygen (and no attached –OH that would indicate an acid).
  5. Return True if any match fully satisfies the acyloin motif, False otherwise.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines whether a molecule (given its SMILES) contains a secondary alpha-hydroxy ketone (acyloin) motif.
    The motif is defined as R–CH(OH)–C(=O)–R′, where:
      - The alpha hydroxy carbon (labeled as :alc) is sp3 with exactly one hydrogen.
      - The carbonyl carbon (labeled as :co) is sp2 and has a double-bonded oxygen and at least one carbon substituent.
      - We also avoid cases where the carbonyl carbon bears a single-bonded oxygen carrying an H (which would indicate an acid).
      
    Args:
      smiles (str): SMILES string representing the molecule.
      
    Returns:
      bool: True if the acyloin motif is found, False otherwise.
      str: Explanation detailing the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern with atom mapping:
    #   [C:alc;H1;X4]         --> the secondary α–hydroxy carbon (with exactly one hydrogen)
    #   ([OX2H])              --> its bonded hydroxyl group (not mapped, we don’t need it later)
    #   -[C:co;X3](=O)[#6]     --> the carbonyl carbon (with a double-bonded oxygen via branch) linked to a carbon substituent.
    acyloin_smarts = "[C:alc;H1;X4]([OX2H])-[C:co;X3](=O)[#6]"
    query = Chem.MolFromSmarts(acyloin_smarts)
    if query is None:
        return False, "Error creating SMARTS pattern"
    
    # Find all substructure matches. Each match is a tuple of 5 atom indices.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Does not contain the required secondary alpha-hydroxy ketone motif"
    
    # Check each match.
    for match in matches:
        # Expecting a tuple of 5 indices; our atoms of interest are:
        #   match[0] -> the alcohol carbon (mapped as :alc)
        #   match[2] -> the carbonyl carbon (mapped as :co)
        if len(match) < 3:
            # Skip if the match is unexpectedly small.
            continue
        
        alc_idx = match[0]
        co_idx = match[2]
        alc_atom = mol.GetAtomWithIdx(alc_idx)
        co_atom  = mol.GetAtomWithIdx(co_idx)
        
        # Confirm the alcohol carbon is secondary (i.e., has exactly one hydrogen).
        if alc_atom.GetTotalNumHs() != 1:
            continue  # Not a secondary carbon
        
        # Check the carbonyl carbon neighbors.
        has_double_bonded_oxygen = False
        has_extra_oh = False
        for nbr in co_atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(co_atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_double_bonded_oxygen = True
                elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # A single-bonded oxygen that carries hydrogen (–OH) would indicate a free acid.
                    if nbr.GetTotalNumHs() >= 1:
                        has_extra_oh = True
        if not has_double_bonded_oxygen:
            continue  # Carbonyl carbon does not have a double-bonded oxygen
        if has_extra_oh:
            continue  # Likely part of a carboxylic acid rather than a ketone
        
        # If this candidate match meets all criteria, we conclude the molecule contains the acyloin motif.
        return True, "Molecule contains the secondary alpha-hydroxy ketone (acyloin) motif"
    
    # If none of the matches qualify, return False.
    return False, "Does not contain a valid secondary alpha-hydroxy ketone (acyloin) motif"

# Optional main block to test a few examples.
if __name__ == '__main__':
    test_examples = {
        "(S)-phenylacetylcarbinol": "C=1C=CC=CC1[C@@H](C(C)=O)O",
        "(S)-benzoin": "O[C@H](C(=O)c1ccccc1)c1ccccc1",
        "Non-acyloin example": "CC(=O)OC1=CC=CC=C1",  # An ester, not an acyloin.
    }
    
    for name, smi in test_examples.items():
        result, reason = is_secondary_alpha_hydroxy_ketone(smi)
        print(f"{name}: {result} -- {reason}")