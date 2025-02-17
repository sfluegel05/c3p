"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin cation
Definition:
  Any organic cation that is an aglycon of anthocyanin cation; they are oxygenated derivatives
  of flavylium (2-phenylchromenylium).

Note: The previous error was due to using Chem.GetSSSR(mol) which in some versions of RDKit
      returns a vector-like object. We now use rdMolDescriptors.CalcNumRings to address this issue.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Criteria:
      - The molecule must be valid.
      - It must be an organic cation (contain at least one positively charged atom).
      - It must have at least three rings (using rdMolDescriptors.CalcNumRings).
      - It must contain a flavylium core pattern (a simplified 2-phenylbenzopyrylium substructure
        with a positively charged oxygen).
        
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an anthocyanidin cation, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for at least one positively charged atom.
    if not any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms()):
        return False, "Molecule does not have a positive formal charge; it is not a cation"
    
    # Calculate the number of rings using rdMolDescriptors.CalcNumRings.
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 3:
        return False, f"Molecule has only {num_rings} ring(s); at least 3 rings are required for an anthocyanidin skeleton"
    
    # Define a simplified SMARTS pattern for the flavylium core.
    # This pattern tries to capture a 2-phenylbenzopyrylium system with a positively charged oxygen atom.
    flavylium_smarts = "c1cc2[O+](ccc2)c(c1)-c3ccccc3"
    flavylium = Chem.MolFromSmarts(flavylium_smarts)
    if flavylium is None:
        return False, "Internal error: invalid SMARTS pattern for flavylium core"
    if not mol.HasSubstructMatch(flavylium):
        return False, "Flavylium core pattern not found in the molecule"
    
    # If all checks pass, classify as an anthocyanidin cation.
    return True, "Molecule matches anthocyanidin cation criteria (organic cation with a flavylium core and sufficient rings)"

# Example usage:
if __name__ == "__main__":
    # Two test examples (one simplified malvidin structure and one of the glycosylated examples).
    example_smiles = [
        "COc1cc(cc(OC)c1O)-c1[o+]c2cc(O)cc(O)c2cc1O",
        "OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@@H]1O"
    ]
    for s in example_smiles:
        result, reason = is_anthocyanidin_cation(s)
        print("SMILES:", s)
        print("Result:", result)
        print("Reason:", reason)
        print("")