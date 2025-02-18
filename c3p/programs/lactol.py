"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol 
Definition: Cyclic hemiacetals formed by intramolecular addition of a hydroxy group to a carbonyl group.
A lactol has, at its hemiacetal (anomeric) center, one –OH group and one –OR group; in a cyclic system the –OR is part of a ring.
"""

from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    The heuristic implemented is as follows:
      1. The molecule is parsed and explicit hydrogens are added.
      2. Search for an atom that is a carbon and is part of a ring.
      3. Check if this carbon has exactly two oxygen neighbors.
      4. Among those oxygens, one should be an -OH group (i.e. have at least one attached hydrogen)
         and the other should be an ether (i.e. no attached hydrogen) that is also in a ring.
    This pattern is consistent with a cyclic hemiacetal (lactol) center.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if a lactol center is identified, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that we can detect -OH groups reliably
    mol = Chem.AddHs(mol)
    
    # Loop over atoms to look for a candidate lactol center:
    for atom in mol.GetAtoms():
        # Focus on carbon atoms in a ring (possible hemiacetal centers)
        if atom.GetAtomicNum() == 6 and atom.IsInRing():
            # Get oxygen neighbors
            oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
            if len(oxy_neighbors) == 2:
                # Count the number of oxygens that bear at least one hydrogen (–OH group)
                oh_count = 0
                or_count = 0
                for oxy in oxy_neighbors:
                    # If the oxygen has any explicit hydrogens, count as –OH.
                    # Note: GetTotalNumHs() includes both implicit and explicit H's if hydrogen counts have been computed.
                    if oxy.GetTotalNumHs() > 0:
                        oh_count += 1
                    else:
                        or_count += 1
                        # additionally, for a lactol the ether oxygen should lie in a ring (intramolecular link)
                        if not oxy.IsInRing():
                            # If the oxygen is not in a ring, then this is likely not an intramolecular hemiacetal.
                            or_count = -99
                            break
                if oh_count == 1 and or_count == 1:
                    return True, f"Found cyclic hemiacetal center at carbon atom index {atom.GetIdx()}"
                    
    return False, "No cyclic hemiacetal (lactol) center found"
    
# Example usage (optional):
if __name__ == "__main__":
    # beta-ascarylopyranose is a lactol (anomeric form in sugars)
    smiles_example = "C[C@@H]1O[C@H](O)[C@H](O)C[C@H]1O"
    result, reason = is_lactol(smiles_example)
    print(result, reason)