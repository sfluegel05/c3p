"""
Classifies: CHEBI:28802 flavonols
"""
"""
Classifies: flavonols
Definition: Any hydroxyflavone in which the ring hydrogen at position 3 
of the heterocyclic (chromen-4-one) ring is replaced by a hydroxy group.
(Also, 3-O–glycosides are accepted.)
The approach here is to search for a 3-hydroxyflavone core.
We use a SMARTS query that encodes a 2-phenylchromen-4-one framework
with an explicit mapping ([O:3]) on the heterocycle carbon at position 3,
which by definition must bear a hydroxy (or glycosylated hydroxy) group.
"""

from rdkit import Chem

def is_flavonols(smiles: str):
    """
    Determines if a molecule is a flavonol (3-hydroxyflavone or 3-O–glycoside)
    based on its SMILES string.
    
    The method uses a SMARTS query that encodes the flavonol core:
      c1ccc(-c2oc3c([O:3])ccc3c2=O)cc1
    This query represents a benzene ring attached to a chromen-4-one system 
    where the heterocyclic carbon at position 3 bears a substituent mapped as [O:3].
    
    After matching the substructure, we inspect the mapped oxygen (position 3):
     - If that oxygen is free (only connected to the flavonol core) it is a free OH.
     - If that oxygen has one extra neighbor, and that neighbor is a methyl group (–CH3),
       we assume it is O-methylated (and the molecule is rejected).
     - If that oxygen has extra neighbors (for instance, in glycosides) then it is accepted.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a flavonol, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS query for 3-hydroxyflavone (flavonol) core.
    # This query encodes a benzene ring (c1ccc…cc1) attached to a chromen-4-one
    # where the heterocyclic ring has an –OH (or O-substituent) at the 3-position.
    query_smarts = "c1ccc(-c2oc3c([O:3])ccc3c2=O)cc1"
    query = Chem.MolFromSmarts(query_smarts)
    if query is None:
        return False, "Error generating SMARTS query for flavonol core"
    
    # Find all substructure matches for the query in the molecule.
    matches = mol.GetSubstructMatches(query, useChirality=True)
    if not matches:
        return False, "Molecule does not contain a valid 3-hydroxyflavone (flavonol) core"
    
    # Map the atom in the query that has the label [O:3].
    mapnum_to_query_idx = {}
    for atom in query.GetAtoms():
        amap = atom.GetAtomMapNum()
        if amap:
            mapnum_to_query_idx[amap] = atom.GetIdx()
    if 3 not in mapnum_to_query_idx:
        return False, "SMARTS query did not map the 3-OH group properly"
    
    # Use the first matching instance.
    match = matches[0]
    flavonol_oh_query_idx = mapnum_to_query_idx[3]
    flavonol_oh_atom_idx = match[flavonol_oh_query_idx]
    oh_atom = mol.GetAtomWithIdx(flavonol_oh_atom_idx)
    
    # Check that this atom is oxygen
    if oh_atom.GetAtomicNum() != 8:
        return False, "The mapped atom for the 3-OH is not oxygen"
    
    # Identify the neighbor in the flavonol core (should be a carbon atom in the fused ring)
    core_neighbors = []
    for nbr in oh_atom.GetNeighbors():
        # Assume the flavonol core carbon(s) belong to rings (at least one ring membership)
        if nbr.IsInRing():
            core_neighbors.append(nbr)
    if not core_neighbors:
        return False, "3-OH oxygen does not appear connected to the flavonol core"
    
    # Now check extra substituents on the hydroxyl oxygen (beyond the flavonol core)
    # Extra substituents (if any) could define a glycoside: if there is one and it is a CH3, then reject.
    extra_neighbors = [nbr for nbr in oh_atom.GetNeighbors() if nbr.GetIdx() not in [atom.GetIdx() for atom in core_neighbors]]
    
    if len(extra_neighbors) == 0:
        return True, "Molecule contains a valid 3-hydroxyflavone core with a free hydroxyl at position 3"
    elif len(extra_neighbors) == 1:
        nbr = extra_neighbors[0]
        # Check if the extra neighbor is a methyl group (CH3). 
        # Methyl is carbon with exactly 3 bound hydrogens and no other heavy atoms attached (besides the OH oxygen).
        if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() == 3:
            return False, "The 3-O substituent appears to be a methoxy group rather than a hydroxyl/glycoside"
        else:
            return True, "Molecule contains a valid 3-hydroxyflavone core with a substituent at position 3"
    else:
        # More than one extra neighbor is consistent with glycosylation or larger substituents.
        return True, "Molecule contains a valid 3-hydroxyflavone core with a glycoside or complex substituent at position 3"

# Example usage (for testing):
if __name__ == "__main__":
    # Example SMILES for tambulin (should be a flavonol)
    test_smiles = "COc1ccc(cc1)-c1oc3c(OC)c(OC)cc(O)c3c(=O)c1O"
    classification, explanation = is_flavonols(test_smiles)
    print("SMILES:", test_smiles)
    print("Classification:", classification)
    print("Reason:", explanation)