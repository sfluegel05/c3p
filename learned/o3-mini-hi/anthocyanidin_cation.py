"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
"""
Classifies: Anthocyanidin cation
Definition:
  Any organic cation that is an aglycon of anthocyanin cation; they are oxygenated derivatives
  of flavylium (2-phenylchromenylium).

This version relaxes the flavylium core SMARTS pattern to be more tolerant of substituents.
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
      - It must contain a flavylium core pattern --
        a fused bicyclic (benzopyrylium) system with a positively charged oxygen and
        an aromatic (phenyl) ring attached at the 2-position.
        
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
    
    # Check that the molecule is an organic cation: at least one atom has a positive formal charge.
    if not any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms()):
        return False, "Molecule does not have a positive formal charge; it is not a cation"
    
    # Calculate the number of rings.
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    if num_rings < 3:
        return False, f"Molecule has only {num_rings} ring(s); at least 3 rings are required for an anthocyanidin skeleton"
    
    # Define a relaxed SMARTS pattern for the flavylium core.
    # This pattern attempts to capture a 2-phenylchromenylium system:
    # - A fused bicyclic aromatic system (benzopyrylium) with the oxygen in the heterocycle 
    #   bearing a positive charge.
    # - A phenyl group attached (at position “2”).
    # The pattern is designed to be tolerant toward additional substituents.
    flavylium_smarts = "[c;R]1ccc2c(c1)[o+][c]([cH]2)-c3ccccc3"
    flavylium = Chem.MolFromSmarts(flavylium_smarts)
    if flavylium is None:
        return False, "Internal error: invalid SMARTS pattern for flavylium core"
    
    if not mol.HasSubstructMatch(flavylium):
        return False, "Flavylium core pattern not found in the molecule"
    
    # If all checks pass, classify as an anthocyanidin cation.
    return True, "Molecule matches anthocyanidin cation criteria (organic cation with a flavylium core and sufficient rings)"

# Example usage: (These examples cover some variants that are known anthocyanidin-derived structures.)
if __name__ == "__main__":
    # A list of example SMILES strings (from anthocyanin derivatives)
    example_smiles = [
        "COc1cc(cc(OC)c1O)-c1[o+]c2cc(O)cc(O)c2cc1O",  # malvidin (simplified)
        "OC[C@H]1O[C@@H](Oc2cc3c(O)cc(O)cc3[o+]c2-c2ccc(O)c(O)c2)[C@H](O)[C@@H](O)[C@@H]1O"  # a glycosylated delphinidin derivative
    ]
    for s in example_smiles:
        result, reason = is_anthocyanidin_cation(s)
        print("SMILES:", s)
        print("Result:", result)
        print("Reason:", reason)
        print("")