"""
Classifies: CHEBI:17478 aldehyde
"""
#!/usr/bin/env python
"""
Classifies: Aldehyde – a compound containing the functional group RC(=O)H 
(plus formaldehyde HCHO) where the carbonyl carbon is bonded to at least one hydrogen.
The substructure is required to have a double-bonded oxygen and exactly two other substituents:
either two hydrogens (formaldehyde) or one hydrogen and one non-oxygen atom (the R group).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is defined as a compound containing the functional group –CHO
    (i.e. a carbonyl group with at least one hydrogen attached and one R group).
    
    To improve on the previous approach we:
      - Generate an aldehyde match using the SMARTS "[CX3;H1,H2](=O)" with explicit hydrogens.
      - For each match, find the double-bonded oxygen and then inspect the two remaining neighbors.
      - Accept a match if (a) the two other neighbors are two hydrogens (formaldehyde), or 
        (b) one is hydrogen and the other is any atom except oxygen (to avoid -COOH, esters, etc.).
    
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
    
    # Add explicit hydrogens for reliable connectivity inspection.
    mol = Chem.AddHs(mol)
    
    # Define the SMARTS pattern for a carbonyl carbon with one or two hydrogens.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3;H1,H2](=O)")
    if aldehyde_pattern is None:
        return False, "Failed to create aldehyde pattern"
    
    # Get initial substructure matches.
    initial_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not initial_matches:
        return False, "No aldehyde group ([CX3;H1,H2](=O)) found in the molecule"
    
    valid_matches = 0
    # Process each potential match.
    for match in initial_matches:
        carbon = mol.GetAtomWithIdx(match[0])
        neighbors = carbon.GetNeighbors()
        dO = None  # to store the double-bonded oxygen (the carbonyl oxygen)
        non_oxo_neighbors = []
        
        # Identify the double-bonded oxygen and collect the others.
        for nbr in neighbors:
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), nbr.GetIdx())
            if bond is None:
                continue
            # Check for double bond and that neighbor is oxygen.
            if bond.GetBondTypeAsDouble() == 2.0 and nbr.GetAtomicNum() == 8:
                # If more than one oxygen is double-bonded (which would be unusual), skip this match.
                if dO is not None:
                    dO = None
                    break
                dO = nbr
            else:
                non_oxo_neighbors.append(nbr)
        
        # The match must have exactly one double-bonded oxygen.
        if dO is None:
            continue
        
        # For a typical aldehyde the carbonyl carbon (sp2)
        # should have exactly two other neighbors (besides the double-bonded oxygen).
        if len(non_oxo_neighbors) != 2:
            continue
        
        # Determine the composition of the two other neighbors.
        # Count how many are hydrogens and capture the non-hydrogen (R) atom.
        h_count = 0
        r_atom = None
        for nbr in non_oxo_neighbors:
            if nbr.GetAtomicNum() == 1:
                h_count += 1
            else:
                r_atom = nbr
        
        # Check valid connectivity:
        # Case 1: Formaldehyde (HCHO) -> two hydrogens.
        # Case 2: Typical aldehyde: exactly one hydrogen and one R group.
        # However, if the R group is oxygen then this is likely an acid/ester and not a true aldehyde.
        if (h_count == 2) or (h_count == 1 and r_atom is not None and r_atom.GetAtomicNum() != 8):
            valid_matches += 1

    if valid_matches == 0:
        return False, "No valid aldehyde group found after connectivity checks"
    
    reason = f"Found {valid_matches} aldehyde group{'s' if valid_matches > 1 else ''} matching the pattern [CX3;H1,H2](=O) with proper connectivity (one carbonyl O and one R group)"
    return True, reason

# Optional testing code; remove or comment out when using in a larger project.
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
        "Cc1cc(O)c(C=O)c2Oc3cc4OC(C)(C)CC(=O)c4c(C)c3OC(=O)c12",  # mollicellin B
        "C[C@H](CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CCC4=CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C)C=O",  # (25R)-3-ketocholest-4-en-26-al
        "OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC=O",  # 20-oxoarachidonic acid
        "[H]C(=O)c1cc(N)ccc1O",  # 5-aminosalicylaldehyde (again)
        "[H]C(=O)c1cccc(O)c1O",  # 2,3-dihydroxybenzaldehyde (again)
        "C(C#C)=O",  # prop-2-ynal should classify as aldehyde
        "[H]C(=O)",  # formaldehyde with explicit hydrogen
        "CC(C)CC=O",  # 3-methylbutanal
        "[H]C(C)=NN(C)C([H])=O",  # gyromitrin (should be accepted)
        "O=C1NC(=O)C2=CC=CC=C12",  # Example that should not be an aldehyde (cyclic imide)
        "O=C1O[C@H](OC[C@@]2(CC3=C([C@@]4(O)CC[C@]4(C[C@@H]3C2)C)C)C)[C@]56[C@]1(C[C@@H]7CC(C)(C)CC7=C5C=O)C6"  # Hexagonin D, a false positive in previous run
    ]
    
    for sm in test_smiles:
        result, explanation = is_aldehyde(sm)
        print(f"SMILES: {sm}\nResult: {result}\nExplanation: {explanation}\n")