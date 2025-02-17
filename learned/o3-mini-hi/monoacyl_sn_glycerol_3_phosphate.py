"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: monoacyl-sn-glycerol 3-phosphate

Definition:
  An sn-glycero-3-phosphate compound having a single unspecified acyl group
  at either position 1 or position 2. In a free glycerol-3-phosphate,
  all three hydroxyl groups are present. In the monoacyl case, one (and only one)
  of the sn-1/sn-2 hydroxyl groups is replaced by an acyl ester group, while the 
  third position (sn-3) is esterified by the phosphate group. 
  
This code uses two explicit SMARTS patterns – one for acylation at sn-1 and 
one for acylation at sn-2 – and accepts the molecule if exactly one of them matches.
It also checks that the molecule contains exactly one phosphorus.
"""

from rdkit import Chem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if the molecule given by the SMILES string is a monoacyl-sn-glycerol 3-phosphate.
    
    A monoacyl-sn-glycerol 3-phosphate is defined as an sn-glycero-3-phosphate compound having
    exactly one acyl group (esterified fatty acid) at either position 1 or 2. In addition,
    the phosphate (attached to position 3) must be monoesterified (i.e. as P(=O)(O)(O)).
    
    The method uses two SMARTS queries for the glycerol-phosphate backbone:
      • For sn-1 acylation: "C(OC(=O)*)C(O)COP(=O)(O)O"
      • For sn-2 acylation: "C(O)C(OC(=O)*)COP(=O)(O)O"
    Only one of these substructure matches must be found.
    
    Also, it requires that there is exactly one phosphorus atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the criteria.
        str: Explanation for the classification decision.
    """
    # 1. Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 2. Ensure exactly one phosphorus is present.
    phosphorus_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if len(phosphorus_atoms) != 1:
        return False, f"Molecule must contain exactly one phosphorus atom (found {len(phosphorus_atoms)})"
    
    # 3. Define two SMARTS patterns corresponding to the glycerol-3-phosphate backbone
    #    with acylation at either sn-1 or sn-2. The "*" (wildcard) in the acyl moiety
    #    makes the chain unspecified.
    #    Pattern for sn-1 acylation: the first carbon of glycerol bears the acyl group.
    sn1_smarts = "C(OC(=O)*)C(O)COP(=O)(O)O"
    #    Pattern for sn-2 acylation: the second carbon of glycerol bears the acyl group.
    sn2_smarts = "C(O)C(OC(=O)*)COP(=O)(O)O"
    
    pattern_sn1 = Chem.MolFromSmarts(sn1_smarts)
    pattern_sn2 = Chem.MolFromSmarts(sn2_smarts)
    
    # 4. Check for matches of each pattern.
    sn1_matches = mol.GetSubstructMatches(pattern_sn1)
    sn2_matches = mol.GetSubstructMatches(pattern_sn2)
    
    # Count the total number of backbone matches.
    total_matches = len(sn1_matches) + len(sn2_matches)
    
    # 5. Accept only if there is exactly one match (i.e. exactly one acylation at sn-1 or sn-2).
    if total_matches == 0:
        return False, "No glycerol-3-phosphate backbone with an acyl substitution (sn-1 or sn-2) found"
    if total_matches > 1:
        return False, f"Multiple ({total_matches}) glycerol-3-phosphate acyl backbone matches found; expected exactly one"
    
    # 6. If we got here, one unique match was found.
    return True, "Molecule has a glycerol-3-phosphate backbone with exactly one acyl ester group at sn-1 or sn-2"

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test a positive case:
    test_smiles = "CCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O"  # 1-icosanoyl-sn-glycero-3-phosphate
    result, reason = is_monoacyl_sn_glycerol_3_phosphate(test_smiles)
    print(result, reason)