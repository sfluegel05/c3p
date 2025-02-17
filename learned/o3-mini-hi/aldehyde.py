"""
Classifies: CHEBI:17478 aldehyde
"""
"""
Classifies: Aldehyde – a compound having the functional group RC(=O)H 
(where the carbonyl carbon is bonded to one hydrogen and one R group or, in the case of formaldehyde, two hydrogens).
The carbonyl group (C=O) must be exocyclic (i.e. not part of a ring) to avoid misclassifying lactones/esters.
"""

from rdkit import Chem

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    
    An aldehyde is defined as a compound containing a carbonyl group (C=O)
    in which the carbonyl carbon is bonded to exactly one hydrogen and one R group 
    (or two hydrogens for formaldehyde). To avoid false positives from cyclic or lactone-type carbonyls,
    we require that the carbonyl double bond is exocyclic.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if at least one valid aldehyde group is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are clear.
    mol = Chem.AddHs(mol)
    
    valid_aldehyde_count = 0
    
    # Iterate over all atoms.
    for atom in mol.GetAtoms():
        # We are only interested in carbon atoms (atomic number 6).
        if atom.GetAtomicNum() != 6:
            continue

        # Identify if this carbon has a double bond to an oxygen,
        # and that C=O bond is exocyclic (i.e. the bond is not in a ring).
        o_dblbond = None  # will store the oxygen atom that is double bonded to this carbon.
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.IsInRing():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:  # oxygen
                    o_dblbond = nbr
                    break
        # If we did not find a qualifying double bonded oxygen, skip this atom.
        if o_dblbond is None:
            continue
        
        # Determine connectivity around this candidate carbon.
        # We use the explicit neighbors (which now include added hydrogens).
        neighbors = atom.GetNeighbors()
        
        # Count hydrogens attached (explicit only).
        h_count = sum(1 for nbr in neighbors if nbr.GetAtomicNum() == 1)
        # Count heavy atoms (atomic number >1) other than the oxygen (already used as the carbonyl partner).
        heavy_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != o_dblbond.GetIdx()]
        
        # Case 1: Typical aldehyde (RC(=O)H): the carbon should have exactly one hydrogen and exactly one heavy R-group.
        if h_count == 1 and len(heavy_neighbors) == 1:
            valid_aldehyde_count += 1
            continue
        # Case 2: Formaldehyde (H2C=O): the carbon should have two hydrogens and no heavy neighbors.
        if h_count == 2 and len(heavy_neighbors) == 0:
            valid_aldehyde_count += 1
            continue
        # Else: not matching the expected connectivity
    if valid_aldehyde_count == 0:
        return False, "No valid aldehyde group found (check connectivity: exocyclic C=O must have one H and one R group, or two H's for formaldehyde)"
    
    # Form an explanatory message.
    group_word = "group" if valid_aldehyde_count == 1 else "groups"
    reason = f"Found {valid_aldehyde_count} aldehyde {group_word} with proper connectivity (exocyclic C=O with one H and one R group, or two H's for formaldehyde)"
    return True, reason

# Optional testing code (can be removed or commented out if not needed)
if __name__ == "__main__":
    # Example SMILES strings (including some from the prompt)
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
        "C(C#C)=O",  # prop-2-ynal
        "[H]C(=O)",  # formaldehyde
        "CC(C)CC=O",  # 3-methylbutanal
    ]
    
    for sm in test_smiles:
        result, explanation = is_aldehyde(sm)
        print(f"SMILES: {sm}\nResult: {result}\nExplanation: {explanation}\n")