"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin cation

Definition:
  Any organic cation that is an aglycon of anthocyanin cation; they are oxygenated derivatives
  of flavylium (2-phenylchromenylium).

This version relaxes the rigid SMARTS matching in favor of a more flexible detection of
a flavylium core via:
  - verifying that the molecule is an organic cation,
  - ensuring the molecule has at least three rings (tricyclic core),
  - and checking that at least one oxygen atom with a +1 formal charge is aromatic, belongs to
    a ring, and is connected to at least one aromatic carbon.
  
Many anthocyanidin derivatives are substituted (glycosylated, acylated, etc.) so an exact
substructure match may fail.
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
      - The molecule is expected to contain a flavylium-like core. In anthocyanidins, an oxygen in the
        central pyran ring typically carries a positive charge and is aromatic and fused into the ring system.
        We check if there is at least one oxygen atom having:
          • a positive formal charge,
          • aromaticity,
          • and residing in a ring and bonded to at least one aromatic carbon.
    
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
    
    # Check that the molecule is an organic cation:
    # We require that at least one atom has a positive formal charge.
    if not any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms()):
        return False, "Molecule does not have a positive formal charge; it is not a cation"
    
    # Check that the molecule has at least three rings (typical for anthocyanidin core).
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 3:
        return False, f"Molecule has {num_rings} ring(s); at least 3 rings are required for an anthocyanidin core"
    
    # Look for a flavylium-like center:
    # Scan for an oxygen atom that is (a) formally positive, (b) aromatic, (c) in a ring, and
    # (d) bonded to at least one aromatic carbon.
    flavylium_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8 and atom.GetFormalCharge() == 1:
            # Check if the oxygen atom is aromatic and is in a ring.
            if not (atom.GetIsAromatic() and atom.IsInRing()):
                continue
            # Check its neighbors: we expect at least one aromatic carbon.
            aromatic_c_neighbor = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIsAromatic():
                    aromatic_c_neighbor = True
                    break
            if aromatic_c_neighbor:
                flavylium_found = True
                break
    
    if not flavylium_found:
        return False, "Flavylium-like center not detected: no aromatic oxygen with +1 charge in ring found"
    
    return True, "Molecule matches anthocyanidin cation criteria (organic cation with a flavylium-like core and sufficient rings)"

# Example usage:
if __name__ == "__main__":
    # A few representative examples (the list provided in the task is extensive)
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