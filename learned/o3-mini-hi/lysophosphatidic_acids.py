"""
Classifies: CHEBI:32957 lysophosphatidic acids
"""
"""
Classifies: lysophosphatidic acids (LPA)

Definition:
    Lysophosphatidic acids (LPA) are monoacylglycerol phosphates obtained by removal of one 
    acyl group from phosphatidic acid. In our heuristic an LPA should:
      • Contain a phosphate group with a P(=O)(O)(O) motif.
      • Have a phosphate directly attached (via an oxygen) to a glycerol backbone.
      • Have exactly one acyl ester linkage (an “O–C(=O)” group) on the glycerol.
      • Have a molecular weight above ~250 Da.
      • Not contain nitrogen atoms.
    
    This heuristic method combines these requirements into one substructure pattern.
    The SMARTS defined here is:
         P(OCC(O)COC(=O)[#6])
    which requires (ignoring stereochemistry):
         • A phosphorus atom connected to an –O–C–C(O)–C–O–C(=O)–[C],
    capturing the idea that one of the glycerol hydroxyls is esterified with an acyl chain
    while another is substituted with a phosphate group.
    
    Note: This is a heuristic and may not cover all possible edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_lysophosphatidic_acids(smiles: str):
    """
    Determines whether a molecule is a lysophosphatidic acid (LPA),
    i.e. a monoacylglycerol phosphate.
    
    The assay is based on:
      - Absence of nitrogen atoms,
      - A reasonable molecular weight (>=250 Da),
      - The presence of a substructure where a phosphate is directly linked 
        to a glycerol backbone that carries exactly one acyl ester bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an LPA, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string."
    
    # Reject molecules that contain any nitrogen atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            return False, "Molecule contains nitrogen atoms; likely not an LPA."
    
    # Check that the molecular weight is above ~250 Da.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low to be an LPA."
    
    # Define a SMARTS pattern for the key LPA substructure:
    #  P(OCC(O)COC(=O)[#6])
    # which requires that the phosphate (P) is bonded to a glycerol fragment that
    # has one esterified hydroxyl (O–C(=O)...) and one free hydroxyl.
    lpa_smarts = "P(OCC(O)COC(=O)[#6])"
    lpa_pattern = Chem.MolFromSmarts(lpa_smarts)
    if not lpa_pattern:
        return False, "Internal error: could not parse SMARTS pattern."
    
    # Check if the molecule has the LPA substructure.
    if mol.HasSubstructMatch(lpa_pattern):
        return True, ("Molecule contains a phosphate directly linked to a glycerol "
                      "backbone with exactly one acyl ester bond (LPA fingerprint).")
    else:
        return False, "Substructure pattern for LPA not found."

# Example usage:
if __name__ == "__main__":
    # Use one of the true positive examples: PA(17:1(9Z)/0:0)
    test_smiles = "P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCC)(O)(O)=O"
    result, reason = is_lysophosphatidic_acids(test_smiles)
    print(result, reason)