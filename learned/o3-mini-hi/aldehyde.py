"""
Classifies: CHEBI:17478 aldehyde
"""
#!/usr/bin/env python
"""
Classifies: Aldehyde – a compound containing the functional group RC(=O)H 
where a carbonyl carbon is bonded to at least one hydrogen (including formaldehyde)
as required by the definition of an aldehyde.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is defined as a compound containing the functional group –CHO 
    (i.e. a carbonyl group with at least one hydrogen attached and one R group).
    
    This function not only finds atoms that match the SMARTS pattern 
    "[CX3;H1,H2](=O)" (a carbonyl carbon with 1 or 2 hydrogens) but then also
    verifies that the only neighbor (aside from the oxygen of the carbonyl double bond)
    is a carbon (or a hydrogen if formyl) so that it isn’t part of an ester, carboxylate, etc.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains at least one valid aldehyde substructure, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to the molecule.
    mol = Chem.AddHs(mol)
    
    # Define an improved aldehyde SMARTS pattern.
    # This pattern accepts a trigonal carbon having 1 or 2 hydrogens (so it catches formaldehyde)
    # and double-bonded to an oxygen.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3;H1,H2](=O)")
    if aldehyde_pattern is None:
        return False, "Failed to create aldehyde pattern"
    
    # Get substructure matches from the molecule.
    initial_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not initial_matches:
        return False, "No aldehyde group ([CX3;H1,H2](=O)) found in the molecule"
    
    valid_matches = 0
    # Evaluate each match with extra checks on the carbon connectivity.
    for match in initial_matches:
        # match[0] is the carbon atom in the supposed aldehyde group.
        carbon = mol.GetAtomWithIdx(match[0])
        
        # For the matched carbon, we need to identify the double-bonded oxygen (the carbonyl oxygen)
        # and then inspect the remaining neighbor(s).
        carbon_neighbors = carbon.GetNeighbors()
        dO = None  # the oxygen involved in the double bond
        other_neighbors = []
        for neighbor in carbon_neighbors:
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
            # Check if this bond is a double bond and neighbor is oxygen.
            if bond.GetBondTypeAsDouble() == 2.0 and neighbor.GetAtomicNum() == 8:
                dO = neighbor
            else:
                other_neighbors.append(neighbor)
        
        # We require exactly one double bonded oxygen.
        if dO is None:
            continue  # skip this match if no proper carbonyl oxygen is found
        
        # Now, for a correct aldehyde:
        # (a) if the carbon is formaldehyde then both other neighbors are hydrogens
        # OR (b) one of the other neighbors is hydrogen and the other is a carbon.
        # (Note: the carbon in an sp2 environment should have three neighbors total.)
        h_count = 0
        carbon_count = 0
        for nbr in other_neighbors:
            if nbr.GetAtomicNum() == 1:
                h_count += 1
            elif nbr.GetAtomicNum() == 6:
                carbon_count += 1
            else:
                # If the neighbor is not hydrogen or carbon (for example, oxygen or N),
                # then it is likely part of an ester/amide, so not a typical aldehyde.
                h_count = -100  # flag as invalid
                break
        
        # We expect either:
        # Formaldehyde: 2 hydrogens (h_count == 2 and carbon_count == 0)
        # or a typical aldehyde: 1 hydrogen and 1 carbon (h_count == 1 and carbon_count == 1)
        if (h_count == 2 and carbon_count == 0) or (h_count == 1 and carbon_count == 1):
            valid_matches += 1
    
    if valid_matches == 0:
        return False, "No valid aldehyde group found after connectivity checks"
    
    reason = f"Found {valid_matches} aldehyde group{'s' if valid_matches > 1 else ''} matching the pattern [CX3;H1,H2](=O) with proper connectivity (one carbonyl O and one R group)"
    return True, reason

# (Optional) Testing code; remove or comment out when using in a larger project.
if __name__ == "__main__":
    test_smiles = [
        "O=CC(CCC=C(C)C)C",  # 5-Heptenal, 2,6-dimethyl-
        "Oc1c(C=O)ccc2ccccc12",  # 1-hydroxy-2-naphthaldehyde
        "CC(=O)O[C@H]1CC[C@]2(C=O)[C@H]3CC[C@@]4(C)[C@@H](CCC4=O)[C@@H]3CC=C2C1",  # 3beta-Hydroxy-17-oxoandrost-5-en-19-al acetate
        "C(CCCCCCCCCCC(=O)[H])CCCCC",  # heptadecanal
        "O=CCCCCCCCCCCCCCCCCCCCCCCC",  # tetracosanal
        "CCCC/C=C/C=O",  # (E)-hept-2-enal
        "CCCCCCCCCCCCCCCCC=O",  # octadecanal
        "[H]C(=O)c1cc(N)ccc1O",  # 5-aminosalicylaldehyde
        "[H]C(=O)c1cccc(O)c1O",  # 2,3-dihydroxybenzaldehyde
        "C(C#C)=O",  # prop-2-ynal (should classify as aldehyde)
        "[H]C(=O)",  # formaldehyde (written with explicit hydrogen)
        "CC(C)CC=O",  # 3-methylbutanal
    ]
    
    for sm in test_smiles:
        result, explanation = is_aldehyde(sm)
        print(f"SMILES: {sm}\nResult: {result}\nExplanation: {explanation}\n")