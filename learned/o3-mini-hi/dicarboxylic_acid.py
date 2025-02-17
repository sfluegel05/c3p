"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Any carboxylic acid containing two carboxy groups (dicarboxylic acid)

This improved version:
  - Adds explicit hydrogens so that –OH is visible.
  - Uses two SMARTS patterns (protonated and deprotonated) for acid groups,
    then post-filters each match: the acid carbon must not be directly bound
    to any nitrogen (to avoid matching carboxamide bonds in peptides).
  - Counts only unique matches based on the acid carbon.
  - Only if exactly two valid carboxy groups are found is the molecule considered a dicarboxylic acid.
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is defined as any molecule containing exactly two carboxy groups.
    
    To avoid false positives (e.g. peptides with amide bonds that have acid groups as part of
    their terminal groups) we filter out any acid pattern where the acid carbon (first atom)
    is directly bound to a nitrogen atom (other than the oxygens that define the acid).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a dicarboxylic acid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the -OH atom is visible (important for matching)
    mol = Chem.AddHs(mol)

    # Define SMARTS patterns for carboxylic acid groups:
    # Protonated acid group: -C(=O)OH
    pat_prot = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    # Deprotonated acid group: -C(=O)[O-]
    pat_deprot = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if pat_prot is None or pat_deprot is None:
        return False, "Error in generating SMARTS patterns for carboxylic acid groups"
    
    # Get all substructure matches for both patterns
    prot_matches = mol.GetSubstructMatches(pat_prot)
    deprot_matches = mol.GetSubstructMatches(pat_deprot)
    
    # We will collect unique “acid centers” defined by the carbon atom (the first atom in each match)
    acid_carbons = set()
    
    # Function: check that the acid carbon is not directly bound to any nitrogen (outside the matched O atoms)
    def valid_acid_match(match):
        acid_idx = match[0]  # in our SMARTS, the first atom is the acid carbon
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        # Go over neighbors: ignore oxygens (atomic num 8) since they are part of the acid group.
        for neighbor in acid_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:
                # found a nitrogen bond that is not part of the standard acid group; skip this match
                return False
        return True

    # Process protonated matches
    for match in prot_matches:
        if valid_acid_match(match):
            acid_carbons.add(match[0])
    
    # Process deprotonated matches
    for match in deprot_matches:
        if valid_acid_match(match):
            acid_carbons.add(match[0])
    
    num_acid_groups = len(acid_carbons)

    # Now decide based on the number of valid acid centers
    if num_acid_groups == 2:
        return True, "Molecule contains exactly two carboxylic acid groups"
    elif num_acid_groups < 2:
        return False, f"Found only {num_acid_groups} carboxylic acid group(s); need exactly two."
    else:
        return False, f"Found {num_acid_groups} carboxylic acid groups; dicarboxylic acid must have exactly two."

# Example test cases (for manual testing)
if __name__ == "__main__":
    examples = {
        "fumarylacetic acid": "OC(=O)CC(=O)\\C=C\\C(O)=O",
        "nocardicin A": "N[C@H](CCOc1ccc(cc1)C(=N\\O)\\C(=O)N[C@H]1CN([C@@H](C(O)=O)c2ccc(O)cc2)C1=O)C(O)=O",
        "meso-tartaric acid": "O[C@@H]([C@@H](O)C(O)=O)C(O)=O",
        "3-ethylmalic acid": "CCC(C(O)C(O)=O)C(O)=O",
        "(4Z)-2-oxohept-4-enedioic acid": "C(=C\\CC(C(O)=O)=O)\\CC(O)=O",
        "Ile-Ser-Asp (peptide false positive?)": "O=C(N[C@H](C(=O)N[C@@H](CC(O)=O)C(O)=O)CO)[C@@H](N)[C@H](CC)C",
    }
    for name, smi in examples.items():
        result, reason = is_dicarboxylic_acid(smi)
        print(f"NAME: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*60}")