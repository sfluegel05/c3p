"""
Classifies: CHEBI:17408 monoacylglycerol
"""
"""
Classifies: Monoacylglycerol
Definition: A glyceride in which exactly one of the three –OH groups of a glycerol backbone is esterified
while the other two remain as free hydroxyl groups (or alkylated), i.e. a monoacylglycerol (MAG).
We require that:
  - The molecule shows exactly one ester group overall.
  - There is a three‐carbon, acyclic (non‐ring) backbone that contains two free –OH substituents and one ester –O–C(=O)R.
We check three possible positions (sn1, sn2 and sn3) where the acyl group can be placed.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_monoacylglycerol(smiles: str):
    """
    Determines if a molecule is a monoacylglycerol based on its SMILES string using a set of refined heuristics.
    
    Our approach:
      1. Parse the SMILES to get the molecule.
      2. Count the total number of ester groups using a generic ester pattern.
         (Monoacylglycerols must have exactly one ester group.)
      3. Use three SMARTS patterns that directly target an acyclic three-carbon (glycerol) backbone
         where exactly one hydroxyl is replaced with an acyl moiety.
         We add the R0 qualifier in these patterns so that the backbone carbons are not part of any ring.
         Patterns for acyl substitution at sn-1, sn-2, or sn-3 are provided.
      4. Also enforce a lower bound on molecular weight (e.g. >140 Da) so that very small molecules are excluded.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a monoacylglycerol, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count total number of ester groups using a general ester pattern.
    # This pattern finds oxygen connected to a carbonyl carbon.
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    num_esters = len(ester_matches)
    if num_esters != 1:
        return False, f"Number of ester groups is {num_esters}; expected exactly 1 for a monoacylglycerol"
    
    # Define three SMARTS patterns for an acyclic glycerol backbone with exactly one acyl (ester) substitution.
    # We require three connected carbons (used non-ring: R0) with appropriate substituents.
    # sn1: acyl substituted at the first terminal carbon.
    sn1_smarts = "O=C(O[CH2;R0])-[CH;R0]-[CH2;R0](O)"
    # sn2: acyl substitution at the central carbon.
    sn2_smarts = "O[CH2;R0](O)-[CH;R0](O C(=O))[CH2;R0](O)"
    # We cannot include a space inside a SMARTS so for sn2 we write it carefully:
    sn2_smarts = "O[CH2;R0](O)-[CH;R0](OC(=O))-[CH2;R0](O)"
    # sn3: acyl substituted at the second terminal carbon.
    sn3_smarts = "O[CH2;R0](O)-[CH;R0]-[CH2;R0](OC(=O))"
    
    # Compile SMARTS patterns.
    pat_sn1 = Chem.MolFromSmarts(sn1_smarts)
    pat_sn2 = Chem.MolFromSmarts(sn2_smarts)
    pat_sn3 = Chem.MolFromSmarts(sn3_smarts)
    
    # Check for a positive match of any of the three patterns.
    backbone_match = (mol.HasSubstructMatch(pat_sn1) or 
                        mol.HasSubstructMatch(pat_sn2) or 
                        mol.HasSubstructMatch(pat_sn3))
    if not backbone_match:
        return False, "Glycerol backbone with one acyl substitution pattern not found"
    
    # Optionally check molecular weight to filter out very small molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 140:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is too low for a monoacylglycerol"
    
    return True, "Contains an acyclic glycerol backbone with exactly one esterified hydroxyl and two free hydroxyl groups, consistent with a monoacylglycerol"

# Testing the function (example test)
if __name__ == "__main__":
    test_smiles = "CCCCCCCC(=O)OC[C@@H](O)CO"  # 1-octanoyl-sn-glycerol
    result, reason = is_monoacylglycerol(test_smiles)
    print(f"SMILES: {test_smiles}\nResult: {result}\nReason: {reason}")