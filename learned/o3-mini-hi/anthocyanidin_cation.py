"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin cation

Definition:
  Any organic cation that is an aglycon of anthocyanin cation; they are oxygenated
  derivatives of flavylium (2-phenylchromenylium).

Improvements over the previous version:
  - We use a relaxed SMARTS pattern for the flavylium core by removing ring-specific
    flags and explicit hydrogen count requirements. This is designed to better match the
    core motif in substituted/glycosylated anthocyanidin structures.
  - We still verify that the molecule has at least one positive formal charge and that it
    has at least three rings.
  
The flavylium core pattern (relaxed):
    "c1ccc2c(c1)[o+][c]c(c2)-c3ccccc3"
This is intended to capture the 2-phenylchromenylium motif.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    
    Criteria:
      - The SMILES string must parse to a valid molecule.
      - The molecule must be an organic cation (contain at least one atom with a positive formal charge).
      - The molecule must have at least three rings.
      - The molecule must contain a relaxed flavylium core pattern that captures a 2-phenylchromenylium
        structure (i.e. a benzopyrylium core with a positively charged oxygen and an aromatic substituent).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an anthocyanidin cation, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is an organic cation by verifying at least one atom has a positive formal charge.
    if not any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms()):
        return False, "Molecule does not have a positive formal charge; it is not a cation"
    
    # Check that the molecule has at least three rings (the anthocyanidin core is tricyclic).
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 3:
        return False, f"Molecule has {num_rings} ring(s); at least 3 rings are required for an anthocyanidin core"
    
    # Define a relaxed SMARTS pattern for the flavylium core.
    # The pattern "c1ccc2c(c1)[o+][c]c(c2)-c3ccccc3" is intended to capture the 2-phenylchromenylium scaffold.
    flavylium_smarts = "c1ccc2c(c1)[o+][c]c(c2)-c3ccccc3"
    flavylium = Chem.MolFromSmarts(flavylium_smarts)
    if flavylium is None:
        return False, "Internal error: invalid SMARTS pattern for flavylium core"
    
    # Check if the molecule contains the flavylium core pattern.
    if not mol.HasSubstructMatch(flavylium):
        return False, "Flavylium core pattern not found in the molecule"
    
    # If all checks pass, return that the molecule qualifies as an anthocyanidin cation.
    return True, "Molecule matches anthocyanidin cation criteria (organic cation with a flavylium core and sufficient rings)"

# Example usage:
if __name__ == "__main__":
    # A short list of representative SMILES strings (these may include derivatives that have additional substituents)
    example_smiles = [
        "COc1cc(cc(OC)c1O)-c1[o+]c2cc(O)cc(O)c2cc1O",  # simplified malvidin
        "OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@@H]1O"  # glycosylated delphinidin derivative
    ]
    for s in example_smiles:
        result, reason = is_anthocyanidin_cation(s)
        print("SMILES:", s)
        print("Result:", result)
        print("Reason:", reason)
        print("")