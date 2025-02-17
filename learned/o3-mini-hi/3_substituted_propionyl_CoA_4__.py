"""
Classifies: CHEBI:65111 3-substituted propionyl-CoA(4-)
"""
"""
Classifies: 3-substituted propionyl-CoA(4-)
Definition:
 'An acyl-CoA(4-) oxoanion arising from deprotonation of the phosphate and diphosphate OH groups
  of any 3-substituted propionyl-CoA; major species at pH 7.3.'

This script uses several structural motifs (SMARTS) to check for:
  - A thioester linkage ([CX3](=O)S) joining the acyl group to the CoA part.
  - An adenine substructure in the CoA headgroup.
  - A characteristic pantetheine fragment ([C@H](O)C(C)(C)COP) which indicates the 3-substituted propionyl center.
  - At least four deprotonated oxygen atoms ([O-]) expected for CoA(4-).
"""

from rdkit import Chem

def is_3_substituted_propionyl_CoA_4__(smiles: str):
    """
    Determines if a molecule matches the class of 3-substituted propionyl-CoA(4-) based on its SMILES string.
    
    Expected features:
      - Thioester group: [CX3](=O)SCC (acyl group bound to CoA).
      - Adenine moiety: our replacement SMARTS "n1cnc2cncnc12" should capture the purine ring system.
      - Pantetheine fragment: [C@H](O)C(C)(C)COP indicating the 3‐substituted center.
      - At least four deprotonated oxygens (atoms with formal charge -1) to account for CoA(4-).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): Tuple with True if classified as 3-substituted propionyl-CoA(4-), otherwise False,
                     and a reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for a thioester linkage.
    # Expected occurrence: the acyl group is connected by a carbonyl to a sulfur followed by at least two carbons.
    thioester_smarts = "[CX3](=O)SCC"
    thioester = Chem.MolFromSmarts(thioester_smarts)
    if not mol.HasSubstructMatch(thioester):
        return False, "Thioester linkage not found"

    # 2. Check for the adenine fragment in the CoA headgroup.
    # We use a more permissive SMARTS ("n1cnc2cncnc12") to catch additional substituents.
    adenine_smarts = "n1cnc2cncnc12"
    adenine = Chem.MolFromSmarts(adenine_smarts)
    if not mol.HasSubstructMatch(adenine):
        return False, "Adenine moiety (CoA head group) not found"

    # 3. Check for the characteristic pantetheine fragment.
    # This motif contains the 3–substituted carbon center.
    pantetheine_smarts = "[C@H](O)C(C)(C)COP"
    pantetheine = Chem.MolFromSmarts(pantetheine_smarts)
    if not mol.HasSubstructMatch(pantetheine):
        return False, "Pantetheine fragment not found – molecular scaffold not matching CoA structure"

    # 4. Check that the molecule has at least four deprotonated oxygen atoms ([O-]).
    deprotonated_oxygens = [atom for atom in mol.GetAtoms() 
                             if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == -1]
    if len(deprotonated_oxygens) < 4:
        return False, f"Found {len(deprotonated_oxygens)} deprotonated oxygens; less than required 4 for CoA(4-)"

    # All tests passed: classify as 3-substituted propionyl-CoA(4-)
    return True, "Molecule matches structural criteria for 3-substituted propionyl-CoA(4-)"

# Example usage:
if __name__ == "__main__":
    # Example SMILES (tetracosanoyl-CoA(4-)):
    smiles_example = ("CCCCCCCCCCCCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)"
                      "[C@H](O)C(C)(C)COP([O-])(=O)OP([O-])(=O)OC[C@H]1O[C@H]"
                      "([C@H](O)[C@@H]1OP([O-])([O-])=O)n1cnc2c(N)ncnc12")
    result, reason = is_3_substituted_propionyl_CoA_4__(smiles_example)
    print(result, reason)