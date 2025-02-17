"""
Classifies: CHEBI:17297 UDP-sugar
"""
#!/usr/bin/env python3
"""
Classifies: UDP-sugar
Definition:
  A pyrimidine nucleotide-sugar having UDP as the nucleotide component attached to
  an unspecified sugar via an anomeric diphosphate linkage.
  
This script uses RDKit to check:
  1) That there is a pyrimidine nucleotide substructure. Here we define it as a uracil moiety that is not methylated (i.e. not thymine).
  2) That there is a diphosphate linkage (P-O-P) with exactly two phosphorus atoms.
If both features are found, the molecule is classified as a UDP-sugar.
"""

from rdkit import Chem

def is_UDP_sugar(smiles: str):
    """
    Determines if a molecule is a UDP-sugar based on its SMILES string.
    
    A UDP-sugar must have:
      - A pyrimidine nucleotide component. Here we demand a uracil ring (n1ccc(=O)[nH]c1=O)
        and reject if a thymine ring is found (which would have a methyl substituent, e.g.
        Cc1c(=O)[nH]c(=O)n1).
      - A diphosphate linkage (a P-O-P connectivity) and exactly 2 phosphorus atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is considered a UDP-sugar, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # We now define a SMARTS for an unmethylated uracil ring.
    # This pattern finds a uracil moiety: (pyrimidine with two carbonyls and no extra substituents on the ring).
    # Note that a thymine (5-methyluracil) would have an extra CH3 at one position.
    uracil_smarts = "n1ccc(=O)[nH]c1=O"
    uracil_pattern = Chem.MolFromSmarts(uracil_smarts)
    if not mol.HasSubstructMatch(uracil_pattern):
        return False, "Uracil (nucleotide) substructure not found"
        
    # Now define a SMARTS for thymine (5-methyluracil).
    # This pattern looks for a methyl group attached to the ring in the position typical for thymine.
    thymine_smarts = "Cc1c(=O)[nH]c(=O)n1"
    thymine_pattern = Chem.MolFromSmarts(thymine_smarts)
    if mol.HasSubstructMatch(thymine_pattern):
        return False, "Thymine (5-methyluracil) substructure found; not UDP-sugar"
    
    # Define SMARTS for a diphosphate (P-O-P) connectivity.
    # This pattern looks for a phosphorus atom connected through an oxygen to another phosphorus atom.
    diphosphate_smarts = "[P]-O-[P]"
    diphosphate_pattern = Chem.MolFromSmarts(diphosphate_smarts)
    if not mol.HasSubstructMatch(diphosphate_pattern):
        return False, "Diphosphate linkage (P-O-P) not found"

    # Count the total number of phosphorus atoms in the molecule.
    n_phosphorus = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if n_phosphorus != 2:
        return False, f"Expected 2 phosphorus atoms in a UDP moiety, found {n_phosphorus}"
        
    return True, "Molecule contains an unmethylated uracil substructure and a diphosphate (P-O-P) linkage indicative of a UDP-sugar"


# Example usage (remove or comment out if integrating into a larger codebase)
if __name__ == "__main__":
    # List of some example SMILES for UDP-sugars (only one example is shown here)
    test_smiles = [
        # UDP-alpha-D-mannuronic acid:
        "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@@H]2O)C(O)=O)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O",
        # This SMILES for dTDP-4-acetamido-4,6-dideoxy-alpha-D-galactose (a known false positive) should be rejected:
        "C[C@H]1O[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H](C[C@@H]2O)n2cc(C)c(=O)[nH]c2=O)[C@H](O)[C@@H](O)[C@H]1NC(C)=O",
        # An example of a UTP molecule (3 phosphates) should also be rejected.
        "O[C@@H]1[C@@H](COP(O)(=O)OP(O)(=O)OP(O)(=O)OCC=C)O[C@H]([C@@H]1O)n1ccc(=O)[nH]c1=O",
    ]
    
    for s in test_smiles:
        result, reason = is_UDP_sugar(s)
        print(f"SMILES: {s}\nResult: {result}\nReason: {reason}\n")