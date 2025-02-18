"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: 3-sn-phosphatidyl-L-serine
Definition:
  A 3-sn-glycerophosphoserine compound having acyl substituents at the 1- and 2-hydroxy positions.
  
Heuristic approach:
1. Parse the SMILES and check for key atoms (phosphorus and nitrogen) that must be found in a phosphoserine.
2. Look for a simplified serine fragment (using the SMARTS "C(N)C(=O)O") so that the head-group is likely present.
3. Identify ester carbonyl centers typical of an acyl ester. We do this by iterating over all carbon atoms
   that are carbonyl centers (i.e. have a double-bonded oxygen) and also are bonded to a single-bonded oxygen.
4. For each such carbonyl we then identify the substituent that is a carbon (and hence likely the acyl chain)
   and compute the length of the carbon-only chain (the longest contiguous path of carbon atoms).
5. Count only those acyl substituents with chain length ≥6. We expect exactly two acyl chains (on sn1 and sn2).
6. Return classification True only if the molecule passes the above tests.
"""

from rdkit import Chem

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    The molecule must contain a glycerophosphoserine head-group (indicated by the presence of P,
    a serine fragment, and by the expected connectivity) and exactly two acyl chains (attached via ester bonds)
    that are sufficiently long (at least 6 carbons).

    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a 3-sn-phosphatidyl-L-serine, else False.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for required atoms: phosphorus (P) must be present and nitrogen (N) to indicate serine.
    atoms = list(mol.GetAtoms())
    if not any(atom.GetAtomicNum() == 15 for atom in atoms):
        return False, "Missing phosphorus (P) required for phosphoserine head-group"
    if not any(atom.GetAtomicNum() == 7 for atom in atoms):
        return False, "Missing nitrogen (N) required for serine"
    
    # Check for a serine fragment using a simplified SMARTS: C(N)C(=O)O
    serine_smarts = "C(N)C(=O)O"
    serine_pattern = Chem.MolFromSmarts(serine_smarts)
    if serine_pattern is None:
        return False, "Internal error processing serine SMARTS pattern"
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "Serine fragment (C(N)C(=O)O) not found in molecule"
    
    # Helper: calculate longest contiguous carbon chain starting from a given carbon atom.
    # We use depth-first search (not worrying about cycles beyond visited atoms).
    def longest_carbon_chain(atom_idx, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        # If the atom is not carbon, chain length is 0.
        if atom.GetAtomicNum() != 6:
            return 0
        max_length = 1
        visited.add(atom_idx)
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            if nbr.GetAtomicNum() == 6:
                chain_length = 1 + longest_carbon_chain(nbr_idx, visited.copy())
                if chain_length > max_length:
                    max_length = chain_length
        return max_length

    # Now scan for potential acyl ester bonds. In an acyl ester group (R-C(=O)-O-R'),
    # the carbonyl carbon (R-C(=O)) should be:
    #   - bonded via a double bond to an oxygen,
    #   - bonded via a single bond to an oxygen (the ester oxygen).
    # We then assume the acyl chain is the substituent bonded directly to the carbonyl carbon that is a carbon.
    acyl_chain_ids = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Look for a double-bonded oxygen neighbor
        dobnd_oxygen = None
        single_oxygens = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    dobnd_oxygen = nbr
                elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    single_oxygens.append(nbr)
        # We expect a carbonyl carbon to have one double-bonded oxygen and one single-bonded oxygen.
        if dobnd_oxygen is None or not single_oxygens:
            continue
        # From the single-bonded oxygen(s), the acyl chain is attached directly to the carbonyl carbon.
        # In a simple ester group, the carbonyl carbon already bears the acyl substituent.
        # (It has three neighbors: the double-bonded oxygen, the ester oxygen, and the acyl chain carbon.)
        # So find the neighbor(s) that are not oxygen.
        acyl_candidate = None
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                acyl_candidate = nbr
                break
        if acyl_candidate is None:
            continue
        # Compute chain length starting at the acyl candidate.
        chain_length = longest_carbon_chain(acyl_candidate.GetIdx(), set())
        # Require at least 6 carbons for a valid acyl chain.
        if chain_length >= 6:
            # Use the index of the carbonyl carbon as a unique handle for the acyl group.
            acyl_chain_ids.add(atom.GetIdx())
    
    # We expect exactly 2 acyl chains on the glycerol backbone (sn1 and sn2)
    if len(acyl_chain_ids) != 2:
        return False, f"Expected 2 acyl groups with sufficiently long chains; found {len(acyl_chain_ids)}"

    return True, "Molecule contains a phosphoserine head-group and 2 acyl groups with long carbon chains"

# Example usage (you can remove or comment out the lines below):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@H](N)C(O)=O)OC(=O)CCCCCCCCCCCCCCCCC"
    result, reason = is_3_sn_phosphatidyl_L_serine(test_smiles)
    print(result, reason)