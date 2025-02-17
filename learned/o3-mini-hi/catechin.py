"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: Catechins
Definition: Members of the class of hydroxyflavan that have a flavan-3-ol skeleton and its substituted derivatives.
The detection uses a combination of a relaxed substructure SMARTS pattern for a 2-phenyl-3,4-dihydro-2H-chromene-3-ol core and a check
that (at least) one free hydroxyl group (–OH) is attached to a saturated carbon in a ring.
Note: This rule–based approach is heuristic and may be adjusted further based on extra data.
"""

from rdkit import Chem

def is_catechin(smiles: str):
    """
    Determines if a molecule belongs to the catechin class based on its SMILES string.
    
    The algorithm is as follows:
      1. Parse the SMILES and remove stereochemistry to avoid mismatches.
      2. Look for a relaxed substructure pattern for the flavan-3-ol (catechin) core.
         The SMARTS pattern used is:
            "c1ccc(c(c1))C2C[C;H](O)c3ccccc3O2"
         This pattern represents a 2-phenyl-3,4-dihydro-2H-chromene (the catechin core),
         where the “C[C;H](O)” part is expected to have a free –OH (at the 3–position).
      3. In addition, we iterate over all oxygen atoms in the molecule and require that at least one
         “free hydroxyl” exists – defined here as an oxygen atom (atomic number 8) that:
           • has only one non‐hydrogen neighbor,
           • that neighbor is an sp3 carbon,
           • and that carbon is in a ring (as expected for the catechin C–ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as catechin, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Remove stereochemistry to simplify matching
    Chem.RemoveStereochemistry(mol)
    
    # Define a relaxed SMARTS pattern for the flavan-3-ol core.
    # The pattern intends to capture a 2-phenylchroman system in which one saturated carbon bears an -OH.
    # (The pattern does not enforce all substituent details.)
    core_smarts = "c1ccc(c(c1))C2C[C;H](O)c3ccccc3O2"
    core_query = Chem.MolFromSmarts(core_smarts)
    if core_query is None:
        return False, "Failed to create SMARTS pattern for catechin core"
    
    # Search for the catechin-like core in the molecule.
    if not mol.HasSubstructMatch(core_query):
        return False, "Molecule does not contain the expected flavan-3-ol core"
    
    # Now, check that there is at least one free hydroxyl group attached to a saturated carbon in a ring.
    # We define a 'free hydroxyl' as an oxygen atom that:
    #  - has atomic number 8,
    #  - is not connected to more than one heavy atom (so not part of an ether or ester),
    #  - and is attached to a carbon in a ring that is sp3 (as expected at C3 of the flavan-3-ol core).
    free_oh_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        # Count heavy-atom (non-hydrogen) neighbors.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            continue
        neighbor = heavy_neighbors[0]
        # Check that the neighbor is a carbon and is sp3.
        if neighbor.GetAtomicNum() == 6 and neighbor.GetHybridization().name == "SP3":
            # And to be strict we require that this carbon is in a ring (the C-ring of flavan-3-ol is ring‐embedded)
            if neighbor.IsInRing():
                free_oh_found = True
                break

    if not free_oh_found:
        return False, "Catechin core detected but no free hydroxyl on a saturated ring carbon was found"
    
    return True, "Molecule contains a flavan-3-ol (catechin) core with at least one free hydroxyl group on a saturated ring carbon"

# Uncomment the section below for simple testing:
# if __name__ == "__main__":
#     test_smiles = "O[C@@H]1Cc2c(O)cc(O)cc2O[C@H]1c1ccc(O)c(O)c1"  # (-)-catechin
#     result, reason = is_catechin(test_smiles)
#     print(result, reason)