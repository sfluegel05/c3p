"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
Definition:
 'An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups 
  of any 3-substituted propionyl-CoA; major species at pH 7.3.'

This script uses several structural motifs (SMARTS) to check for a thioester linkage, the CoA headgroup (adenine)
and the characteristic pantetheine fragment (which includes a 3–substituted carbon center) as well as the proper 
number of deprotonated oxygens.
"""

from rdkit import Chem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule matches the class of 3-substituted propionyl-CoA(4-) based on its SMILES string.
    The expected features are:
      - A thioester group [CX3](=O)S (for the acyl-CoA bond)
      - An adenine substructure (the purine moiety of CoA)
      - A pantetheine fragment containing [C@H](O)C(C)(C)COP which reflects the 3-substituted nature
      - At least 4 deprotonated oxygens (formal charge -1) consistent with CoA(4-)
      
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str): Tuple with the result (True if classified as 3-substituted propionyl-CoA(4-), else False)
                     and a reason for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for a thioester linkage.
    # The acyl group connected via a carbonyl to a sulfur (which then attaches to the CoA part).
    thioester_smarts = "[CX3](=O)SCC"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester):
        return False, "Thioester linkage not found"
    
    # 2. Check for the adenine fragment of the CoA moiety.
    # Many acyl-CoA compounds contain an adenine ring. The SMARTS below should pick up the purine pattern.
    adenine_smarts = "c1ncnc2n(c1)cnc2"
    adenine = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine):
        return False, "Adenine moiety (CoA head group) not found"
    
    # 3. Check for the characteristic pantetheine fragment.
    # This motif contains the 3–substituted carbon center ([C@H](O)C(C)(C)COP) that is present in all CoA’s.
    pantetheine_smarts = "[C@H](O)C(C)(C)COP"
    pantetheine = Chem.MolFromSmarts(pantetheine_smarts)
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Pantetheine fragment not found – molecular scaffold not matching CoA structure"
    
    # 4. Check that the molecule has at least four deprotonated oxygens ([O-])
    deprotonated_oxygens = [atom for atom in mol.GetAtoms() 
                             if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1]
    if len(deprotonated_oxygens) < 4:
        return False, f"Found {len(deprotonated_oxygens)} deprotonated oxygens; less than required 4 for CoA(4-)"
    
    # If all tests are passed then we conclude that the molecule matches our class.
    return True, "Molecule matches structural criteria for 3-substituted propionyl-CoA(4-)"

# Example usage:
if __name__ == "__main__":
    # Example SMILES from the user specifications (tetracosanoyl-CoA(4-))
    smiles_example = ("CCCCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)"
                      "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]"
                      "([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    result, reason = is_3_substituted_propionyl_CoA_4__(smiles_example)
    print(result, reason)