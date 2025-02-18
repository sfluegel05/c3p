"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
"""
Classifies: Secondary Alpha-Hydroxy Ketone (Acyloin)

Definition:
  A secondary α‐hydroxy ketone (acyloin) is defined as a molecule that contains the motif
      R–CH(OH)–C(=O)–R′,
  where the “CH(OH)” center is secondary (i.e. it has exactly one hydrogen and two heavy‐atom neighbors)
  and the carbonyl carbon is bonded to a double-bonded oxygen (and at least one carbon substituent, i.e. not an acid).
  
Strategy:
  1. We start with a SMARTS query that looks for an sp3 carbon with exactly one H attached,
     which bears an –OH and is linked by a single bond to an sp2 carbon with a double-bonded oxygen.
     The chosen query is: "[C;H1;X4]([OX2H])-[C;X3](=O)[#6]".
  2. Then, for every match, we verify:
       • The alcohol (α–hydroxy) carbon actually has exactly one hydrogen.
       • The carbonyl carbon has an attached oxygen via a double bond and does not also bear a single–bonded –OH (to avoid free acid cases).
  3. If any match passes these extra checks we return True with a corresponding reason.
  
Note: There remains a challenge to fully separate “true” acyloins from similar substructures in very complex molecules.
"""

from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines whether a molecule (given its SMILES) contains a secondary alpha-hydroxy ketone (acyloin) motif.
    The motif is defined as R–CH(OH)–C(=O)–R′, where the CH(OH) center is secondary (one hydrogen) and
    the carbonyl carbon is part of a ketone (has a double-bonded oxygen) and is not part of a free acid.
    
    Args:
      smiles (str): SMILES string representing the molecule.
    
    Returns:
      bool: True if the acyloin motif is found, False otherwise.
      str: Explanation detailing the decision.
    """
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS pattern.
    # Explanation:
    #   [C;H1;X4]         : an sp3 carbon with exactly one attached hydrogen (the α–hydroxy carbon)
    #   ([OX2H])          : that carbon is bonded to an –OH group.
    #   -[C;X3](=O)[#6]    : via a single bond to a carbon (in sp2, X3) that bears a double-bonded oxygen and is bonded to at least one carbon.
    acyloin_smarts = "[C;H1;X4]([OX2H])-[C;X3](=O)[#6]"
    query = Chem.MolFromSmarts(acyloin_smarts)
    if query is None:
        return False, "Error creating SMARTS pattern"
    
    # Find all substructure matches.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Does not contain the required secondary alpha-hydroxy ketone motif"
    
    # Process each candidate match.
    for match in matches:
        # In our SMARTS we expect 3 atoms:
        #   match[0] : the α–hydroxy (alcohol) carbon (should be sp3 and secondary)
        #   match[1] : the carbonyl carbon (should have a double-bonded oxygen)
        #   match[2] : one substituent attached to the carbonyl carbon (must be carbon)
        alc_idx, co_idx, sub_idx = match
        alc_atom = mol.GetAtomWithIdx(alc_idx)
        co_atom  = mol.GetAtomWithIdx(co_idx)
        
        # Verify that the alcohol carbon has exactly one hydrogen. (Secondary: only one H)
        if alc_atom.GetTotalNumHs() != 1:
            # Not the correct secondary center; skip this match.
            continue
        
        # Now check that the carbonyl carbon (co_atom) has a proper double-bonded oxygen.
        # We expect one neighbor to be an oxygen with a double bond.
        has_dbl_oxygen = False
        # Also, if the carbonyl carbon has a single bonded oxygen with an H (an –OH), it is more like an acid.
        has_extra_OH = False
        for nbr in co_atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(co_atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    has_dbl_oxygen = True
                elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Check if this oxygen bears an H (–OH) indicating an acid moiety.
                    if nbr.GetTotalNumHs() >= 1:
                        has_extra_OH = True
        if not has_dbl_oxygen:
            # The carbonyl carbon does not have a proper C=O bond.
            continue
        if has_extra_OH:
            # Likely a carboxylic acid rather than a ketone.
            continue
        
        # If we reach here, we have found a candidate acyloin motif.
        return True, "Molecule contains the secondary alpha-hydroxy ketone (acyloin) motif"
    
    # After checking all candidate matches, if none pass the refined criteria, return false.
    return False, "Does not contain a valid secondary alpha-hydroxy ketone (acyloin) motif"

# Optional main section to quickly test a few examples.
if __name__ == '__main__':
    examples = {
        "(S)-phenylacetylcarbinol": "C=1C=CC=CC1[C@@H](C(C)=O)O",
        "dTDP-4-dehydro-2,6-dideoxy-D-glucose": "C[C@H]1OC(C[C@@H](O)C1=O)OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H](C[C@@H]1O)n1cc(C)c(=O)[nH]c1=O",
        "(S)-benzoin": "O[C@H](C(=O)c1ccccc1)c1ccccc1",
        "False positive candidate": "[H][C@]12C[C@@]1([H])[C@]1(C)[C@@H](O)C(=O)\\C(=C(\\C)C(=O)OC)[C@@]3([H])C1=C2C[C@@]1([H])[C@@]2(C)[C@]4([H])C[C@]4([H])[C@](O)(C\\C(C)=C\\C(=O)OC)[C@@]2([H])CC2=C(CO)C(=O)O[C@]312",
    }
    
    for name, smi in examples.items():
        result, reason = is_secondary_alpha_hydroxy_ketone(smi)
        print(f"{name}: {result} -- {reason}")