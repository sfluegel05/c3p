"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: Flavonoid (defined as any member of the 'superclass' flavonoids 
whose skeleton is based on 1-benzopyran with an aryl substituent at position 2).
This function uses a substructure search to identify a 2-phenylbenzopyran
pattern, which is a common core in many flavonoids.
"""

from rdkit import Chem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a (putative) flavonoid based on its SMILES string.
    The classifier searches for a 2-phenylbenzopyran core (a benzopyran motif
    with an aryl substituent at position 2). Note that many flavonoids are decorated
    with additional substituents (e.g., hydroxyl groups, glycosides, methyl ethers)
    and slight variations exist. This implementation uses one heuristic SMARTS query.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule appears to contain a flavonoid core, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the 2-phenylbenzopyran core.
    # Explanation of the SMARTS:
    #   - "c1ccc2c(c1)" matches an aromatic ring fused to a second ring.
    #   - "OC(C2)" represents a pyran (six-membered oxygen heterocycle) where the
    #     carbon (in position “2”) is not aromatic (thus allowing for sp3 stereochemistry, as in flavanones or flavanols).
    #   - "c3ccccc3" requires that a phenyl ring is directly attached to that carbon.
    #
    # In many flavonoids the core is a 2-phenylchroman (or a derivative thereof).
    pattern_smarts = "c1ccc2c(c1)OC(C2)c3ccccc3"
    core_pattern = Chem.MolFromSmarts(pattern_smarts)
    if core_pattern is None:
        return None, None  # If for some reason the pattern is not created

    # Check whether the molecule contains the flavonoid core
    if mol.HasSubstructMatch(core_pattern):
        return True, "Molecule contains a 2-phenylbenzopyran core typical for flavonoids."
    
    # Otherwise, report that the key core was not found
    return False, "2-phenylbenzopyran core not found – molecule does not appear to be a flavonoid."

# Example usage:
if __name__ == "__main__":
    # Testing with (R)-naringenin (a common flavonoid)
    test_smiles = "Oc1ccc(cc1)[C@H]1CC(=O)c2c(O)cc(O)cc2O1"
    result, reason = is_flavonoid(test_smiles)
    print(result, ":", reason)