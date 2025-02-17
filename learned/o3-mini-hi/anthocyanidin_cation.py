"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin cation

Definition:
  Any organic cation that is an aglycon of anthocyanin cation; they are oxygenated derivatives
  of flavylium (2-phenylchromenylium).

Improved criteria:
  - The SMILES string must parse to a valid molecule.
  - The molecule (or major fragment) must have an overall positive formal charge and must not contain
    any negatively charged atoms.
  - The molecule must have at least three rings.
  - The molecule must contain a flavylium-like core, which we detect by a SMARTS match for a fused
    tricyclic 2-phenylchromenylium motif.
  
Note: Because many anthocyanidin derivatives are substituted (glycosylated, acylated, etc.), an exact
SMARTS might fail sometimes. Nevertheless, this version tries to be stricter by excluding molecules that
contain any negative charges and by requiring a match to a flavylium-like SMARTS pattern.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    
    Improved Criteria:
      - The SMILES must produce a valid molecule.
      - The molecule must have an overall positive net formal charge AND have no atoms with a negative charge.
      - The molecule must contain at least three rings.
      - The molecule must contain a fused tricyclic flavylium core.
        We use a SMARTS pattern for a simplified flavylium core:
          "c1ccc2c(c1)[c]([O+])c3ccccc23"
        which looks for a benzopyrylium fused with another benzene ring (representing 2-phenylchromenylium).
    
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
    
    # Compute overall net formal charge and check for any negatives.
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge <= 0:
        return False, f"Overall formal charge is {total_charge}; expected a positive charge for an anthocyanidin cation"
    if any(atom.GetFormalCharge() < 0 for atom in mol.GetAtoms()):
        return False, "Molecule contains one or more atoms with negative formal charge; not a pure organic cation"

    # Check that the molecule has at least three rings (typical anthocyanidin core has 3 fused rings)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 3:
        return False, f"Molecule has {num_rings} ring(s); at least 3 rings are required for an anthocyanidin core"
    
    # To capture a flavylium-like core reliably, use a SMARTS that looks for a fused three-ring motif.
    # This SMARTS tries to detect a benzopyrylium structure: a benzene ring fused to a pyran ring
    # (with an oxygen carrying a + charge) and connected to another aromatic ring.
    flavylium_smarts = "c1ccc2c(c1)[c]([O+])c3ccccc23"
    flavylium_core = Chem.MolFromSmarts(flavylium_smarts)
    if flavylium_core is None:
        return False, "Error creating flavylium SMARTS pattern"
    
    if not mol.HasSubstructMatch(flavylium_core):
        return False, "Flavylium-like core not detected (no match to expected tricyclic motif with an oxygen+)"
    
    return True, "Molecule matches anthocyanidin cation criteria (organic cation with flavylium-like fused core and sufficient rings)"

# Example usage:
if __name__ == "__main__":
    # A few test examples
    test_smiles = [
        # True positive examples (simplified or aglycon-like cores)
        "COc1cc(cc(OC)c1O)-c1[o+]c2cc(O)cc(O)c2cc1O",  # simplified malvidin like structure
        "OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@@H]1O",  # glycosylated delphinidin derivative
        
        # False positive candidates from previous attempt would have negative charges;
        # for example, an anthocyanidin "olate" would be rejected.
        "[Cl-].OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)c3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@@H]1O",
    ]
    
    for s in test_smiles:
        result, reason = is_anthocyanidin_cation(s)
        print("SMILES:", s)
        print("Result:", result)
        print("Reason:", reason)
        print("")