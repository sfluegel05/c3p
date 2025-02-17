"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: Ceramide (N-acyl-sphingoid bases)

Ceramides are defined as sphingoid base derivatives that contain an amide-linked fatty acid.
The fatty acid is typically saturated or monounsaturated with a chain length from 14 to 26 carbons.
A sphingoid base generally has an amino group (whose acylation gives the amide) and at least one hydroxyl group (often at carbon 2).
This implementation uses simple SMARTS patterns and a routine to measure the length of the fatty acyl chain.
Note: This is a heuristic classifier and may not capture every nuance of ceramide structure.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide, based on its SMILES string.
    
    Checks for:
      1. The presence of an amide group (C(=O)N).
      2. That one side of an amide bond carries a contiguous alkyl chain (only carbon atoms)
         having between 14 and 26 atoms (the typical fatty acid chain).
      3. A rough sphingoid base pattern, here defined as a nitrogen attached to a carbon
         that is also bound to a hydroxymethyl group (i.e. pattern "N[C](CO)").
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a ceramide, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. Check for an amide functional group ---
    # Using a simple SMARTS for an amide: carbonyl carbon attached to a nitrogen.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide group (C(=O)N) found"
    
    # --- 2. Check for a sphingoid base pattern ---
    # This is a heuristic: we look for a pattern where a nitrogen is attached to a carbon
    # that in turn is attached to a –CH2OH group. (Many sphingoid bases have an N-C(CO) motif.)
    sphingo_smarts = Chem.MolFromSmarts("N[C](CO)")
    if not mol.HasSubstructMatch(sphingo_smarts):
        return False, "No sphingoid base pattern (e.g. N[C](CO)) found"
    
    # --- Helper: Compute the longest contiguous carbon chain from a starting carbon.
    def longest_chain(mol, current_idx, coming_from_idx, visited):
        """
        Recursively compute the length of the longest contiguous chain 
        (only following carbon atoms) starting from current_idx.
        
        Args:
            mol: RDKit Mol object.
            current_idx: int, current atom index.
            coming_from_idx: int, the atom index from which we came.
            visited: set, indices already visited in the current chain.
            
        Returns:
            int: longest chain length (counting the current carbon as 1).
        """
        max_length = 1
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nb in current_atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx == coming_from_idx:
                continue
            if nb.GetAtomicNum() != 6:  # Follow only carbon atoms
                continue
            if nb_idx in visited:
                continue
            # Copy visited for each branch to avoid cross‐branch contamination.
            new_visited = visited.copy()
            new_visited.add(nb_idx)
            branch_length = 1 + longest_chain(mol, nb_idx, current_idx, new_visited)
            if branch_length > max_length:
                max_length = branch_length
        return max_length

    # --- 3. Check that one amide bond has a fatty acyl chain in the correct range ---
    fatty_acyl_found = False
    for match in amide_matches:
        # In the SMARTS "C(=O)N", match[0] is the carbonyl carbon and match[1] is the nitrogen.
        carbonyl_idx = match[0]
        nitrogen_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        # Look at neighbors of the carbonyl carbon.
        # One neighbor is the nitrogen (already in the amide) and one is the carbonyl oxygen.
        # We look for a carbon neighbor that could be the start of a fatty acid chain.
        fatty_neighbors = []
        for nb in carbonyl_atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx == nitrogen_idx:
                continue
            if nb.GetAtomicNum() == 8:  # likely the carbonyl oxygen
                continue
            if nb.GetAtomicNum() == 6:
                fatty_neighbors.append(nb_idx)
        # For each candidate start, check its contiguous carbon chain length.
        for start_idx in fatty_neighbors:
            chain_length = longest_chain(mol, start_idx, carbonyl_idx, {start_idx})
            # The fatty acid chain length in ceramides is typically 14–26 carbons.
            if 14 <= chain_length <= 26:
                fatty_acyl_found = True
                break
        if fatty_acyl_found:
            break
    if not fatty_acyl_found:
        return False, "No fatty acyl chain (14-26 contiguous carbons) found attached via amide linkage"
    
    return True, "Molecule contains an amide-linked fatty acid (14-26 carbons) and a sphingoid base backbone"