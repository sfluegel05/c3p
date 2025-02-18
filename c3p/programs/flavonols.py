"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: flavonols
Definition: Any hydroxyflavone in which the ring hydrogen at position 3 
of the heterocyclic (chromen-4-one) ring is replaced by a hydroxy group.
(Also, 3-O–glycosides are accepted.)
Examples include tambulin, pinoquercetin, etc.
"""

from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol (3-hydroxyflavone)
    based on its SMILES string.
    
    The approach is to search for a 2-phenylchromen-4-one core (flavone core)
    with a hydroxyl substituent at position 3. To allow for additional substituents
    (such as glycosides, extra methoxies, etc.) on the aromatic rings we use aromatic
    atom queries ("a") in the SMARTS rather than rigid aliphatic ("c").
    
    The SMARTS we use is:
      a1ccc(-a2oc3cc([O:3])cc(=O)c3c2)cc1
    which looks for a benzene ring (ring "1") attached to a heterocycle ("2" and "3")
    where the heterocycle carries an –OH (atom mapped “3”) at the position corresponding
    to C3 in the flavone nomenclature.
    
    After matching the core, the substituent attached to the oxygen mapped as 3
    is inspected. If the oxygen bears no extra neighbor (besides the core carbon), it is 
    assumed to be a free –OH. If it has one additional neighbor that is a CH3 (i.e. a methoxy),
    the molecule is rejected. Otherwise (e.g. a glycoside), the match is accepted.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a flavonol, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS query for the flavonol core.
    # Using aromatic atom query "a" (instead of strict "c") allows extra substituents.
    # The SMARTS enforces a 2-phenylchromen-4-one core with an -OH at C3 (mapped as [O:3]).
    query_smarts = 'a1ccc(-a2oc3cc([O:3])cc(=O)c3c2)cc1'
    query = Chem.MolFromSmarts(query_smarts)
    if query is None:
        return False, "Error in generating SMARTS query for flavonol core"
    
    # Find all substructure matches
    matches = mol.GetSubstructMatches(query, useChirality=True)
    if not matches:
        return False, "Molecule does not contain a valid 3-hydroxyflavone (flavonol) core"
    
    # For interpretation of the match, we need to know which query atom is mapped as 3.
    mapnum_to_query_idx = {}
    for atom in query.GetAtoms():
        amap = atom.GetAtomMapNum()
        if amap:
            mapnum_to_query_idx[amap] = atom.GetIdx()
    if 3 not in mapnum_to_query_idx:
        return False, "SMARTS query did not map the 3-OH group properly"
    
    # Use the first match (if multiple exist, one valid match is enough).
    match = matches[0]
    flavonol_oh_query_idx = mapnum_to_query_idx[3]
    flavonol_oh_idx = match[flavonol_oh_query_idx]
    oh_atom = mol.GetAtomWithIdx(flavonol_oh_idx)
    
    # Identify the core carbon to which the hydroxyl is attached.
    # In our query the OH is attached to exactly one core atom.
    core_neighbor = None
    nbrs = oh_atom.GetNeighbors()
    if not nbrs:
        return False, "3-OH oxygen not attached to any atom (unexpected)"
    else:
        # It should at least be attached to the flavonol carbon in the core
        core_neighbor = nbrs[0]
    
    # Examine extra neighbors (besides the core connection)
    extra_neighbors = [nbr for nbr in oh_atom.GetNeighbors() if nbr.GetIdx() != core_neighbor.GetIdx()]
    
    # Decide based on substituent on the oxygen:
    if len(extra_neighbors) == 0:
        # Free hydroxyl
        return True, "Molecule contains a valid 3-hydroxyflavone core with a free hydroxyl at position 3"
    elif len(extra_neighbors) == 1:
        nbr = extra_neighbors[0]
        # If the extra neighbor is a methyl group (CH3), it indicates a methoxy group.
        if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() == 3:
            return False, "The 3-O substituent appears to be a methoxy group rather than a hydroxyl or glycoside"
        else:
            return True, "Molecule contains a valid 3-hydroxyflavone core with a substituent at position 3"
    else:
        # More than one extra neighbor (as in glycosidic attachments) is acceptable.
        return True, "Molecule contains a valid 3-hydroxyflavone core with a glycoside or larger substituent at position 3"


# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test example: tambulin should be classified as a flavonol.
    test_smiles = "COc1ccc(cc1)-c1oc2c(OC)c(OC)cc(O)c2c(=O)c1O"
    classification, explanation = is_flavonols(test_smiles)
    print("SMILES:", test_smiles)
    print("Classification:", classification)
    print("Reason:", explanation)