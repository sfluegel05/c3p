"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: lysophosphatidic acids (LPA)

Definition:
  Lysophosphatidic acids (LPA) are monoacylglycerol phosphates obtained by removal of one 
  acyl group from phosphatidic acid. Here our heuristic requires:
    • No nitrogen atoms present.
    • A molecular weight roughly between 250 and 600 Da.
    • Exactly one phosphorus atom.
    • A phosphate group (P(=O)(O)(O)) directly linked (via an oxygen) to a glycerol backbone
      that carries exactly one acyl ester linkage.
  Because the acyl group may be attached at different positions on the glycerol backbone, 
  we test against two SMARTS patterns:
    Pattern 1: "P(OCC(O)COC(=O)[#6])"  (acyl group at terminal glycerol carbon)
    Pattern 2: "P(OCC(OC(=O)[#6])CO)"  (acyl group at the middle glycerol position)
  
Note: This heuristic is imperfect but considerably improves classification.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines whether a molecule is a lysophosphatidic acid (LPA),
    i.e. a monoacylglycerol phosphate.
    
    The assay checks that:
      - The molecule does not contain any nitrogen atoms.
      - Its molecular weight is in the range [250,600] Da.
      - It contains exactly one phosphorus atom.
      - It contains one of the two LPA substructure fingerprints.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an LPA, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Check: No nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            return False, "Molecule contains nitrogen atoms; likely not an LPA."
    
    # Check: Molecular weight between 250 and 600 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low to be an LPA."
    if mol_wt > 600:
        return False, "Molecular weight too high; likely not an LPA."
    
    # Check: Exactly one phosphorus atom.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return False, "Molecule must contain exactly one phosphorus atom."
    
    # Define two SMARTS patterns for the LPA fingerprint:
    # Pattern 1: phosphate attached to glycerol with acyl group on the terminal carbon.
    lpa_smarts1 = "P(OCC(O)COC(=O)[#6])"
    # Pattern 2: alternative connectivity (acyl group on the middle glycerol carbon)
    lpa_smarts2 = "P(OCC(OC(=O)[#6])CO)"
    
    pattern1 = Chem.MolFromSmarts(lpa_smarts1)
    pattern2 = Chem.MolFromSmarts(lpa_smarts2)
    if pattern1 is None or pattern2 is None:
        return False, "Internal error: could not parse SMARTS patterns."
    
    # Check for substructure matches.
    match1 = mol.HasSubstructMatch(pattern1)
    match2 = mol.HasSubstructMatch(pattern2)
    
    if match1 ^ match2:  # Exactly one (exclusive OR) of these patterns should be present.
        return True, ("Molecule contains a phosphate directly linked to a glycerol backbone with "
                      "exactly one acyl ester bond (LPA fingerprint).")
    elif match1 and match2:
        # If both patterns match, this is ambiguous.
        return False, "Ambiguous substructure: multiple LPA fingerprints detected."
    else:
        return False, "Substructure pattern for LPA not found."

# Example usage:
if __name__ == "__main__":
    # Test using one of the true positive examples: PA(17:1(9Z)/0:0)
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(O)(O)=O"
    result, reason = is_lysophosphatidic_acids(test_smiles)
    print(result, reason)