"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
"""
Classifies: Any carboxylic acid containing two carboxy groups (dicarboxylic acid)

Improvement notes:
  - Instead of using one SMARTS pattern for –COOH groups,
    we now use two SMARTS patterns:
      • one for the protonated acid: [CX3](=O)[OX2H1]
      • one for the deprotonated acid: [CX3](=O)[O-]
  - We count unique acid “centers” based on the carbon atom in each match.
  - If exactly two distinct acid carbons are found, we classify as a dicarboxylic acid.
  
  Note: This approach may still give unexpected results for very complex molecules,
  especially for peptides or macromolecules that happen to contain exactly two carboxy groups.
"""

from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid based on its SMILES string.
    A dicarboxylic acid is defined as any molecule containing exactly two carboxy groups –
    as recognized by either a protonated ([CX3](=O)[OX2H1]) or deprotonated ([CX3](=O)[O-]) pattern.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a dicarboxylic acid, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for carboxylic acid groups:
    # Pattern for the protonated carboxyl group –COOH
    pat_prot = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    # Pattern for the deprotonated carboxylate –COO–
    pat_deprot = Chem.MolFromSmarts("[CX3](=O)[O-]")
    
    if pat_prot is None or pat_deprot is None:
        return False, "Error in generating SMARTS patterns for carboxylic acid groups"
    
    # Find all substructure matches for each pattern
    prot_matches = mol.GetSubstructMatches(pat_prot)
    deprot_matches = mol.GetSubstructMatches(pat_deprot)
    
    # We will make a set of unique acid centers (the carbon atom index of the [CX3] in the pattern)
    acid_carbons = set()
    for match in prot_matches:
        # In our pattern, the first atom is the acid carbon
        acid_carbons.add(match[0])
    for match in deprot_matches:
        acid_carbons.add(match[0])
    
    num_acid_groups = len(acid_carbons)
    
    # Decision: exactly two acid centers classify as dicarboxylic acid
    if num_acid_groups == 2:
        return True, "Molecule contains exactly two carboxylic acid groups"
    elif num_acid_groups < 2:
        return False, f"Found only {num_acid_groups} carboxylic acid group(s); need exactly two."
    else:
        return False, f"Found {num_acid_groups} carboxylic acid groups; dicarboxylic acid must have exactly two."

# Example test cases
if __name__ == "__main__":
    # A few examples: 
    examples = {
        "fumarylacetic acid": "OC(=O)CC(=O)\\C=C\\C(O)=O",
        "meso-tartaric acid": "O[C@@H]([C@@H](O)C(O)=O)C(O)=O",
        "3-ethylmalic acid": "CCC(C(O)C(O)=O)C(O)=O",
        "(4Z)-2-oxohept-4-enedioic acid": "C(=C\\CC(C(O)=O)=O)\\CC(O)=O",
        "Ile-Ser-Asp (peptide false positive?)": "O=C(N[C@H](C(=O)N[C@@H](CC(O)=O)C(O)=O)CO)[C@@H](N)[C@H](CC)C"
    }
    
    for name, smi in examples.items():
        result, reason = is_dicarboxylic_acid(smi)
        print(f"NAME: {name}\nSMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*60}")