"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
"""
Classifies: Tetrahydrofuranone 
Definition: Any oxolane (i.e. tetrahydrofuran ring, a five‐membered ring 
with one oxygen and four carbons) having an oxo- substituent (exocyclic carbonyl group)
at any position on the tetrahydrofuran ring.
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    The molecule must include at least one 5-membered ring with exactly one oxygen 
    (the oxolane ring) in which one of the ring carbons has an exocyclic double bonded oxygen.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is a tetrahydrofuranone, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # Iterate through all rings in the molecule.
    for ring in atom_rings:
        # Check if ring is 5-membered.
        if len(ring) != 5:
            continue
        
        # Count oxygens and carbons in the ring.
        oxygen_count = 0
        carbon_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxygen_count += 1
            elif atom.GetAtomicNum() == 6:
                carbon_count += 1
            else:
                # If any atom is not C or O, then skip this ring.
                oxygen_count = -100
                break
        # We want exactly one oxygen (oxolane) with four carbons.
        if oxygen_count != 1 or carbon_count != 4:
            continue
        
        # Now check the candidate ring: find a carbon within the ring that has an exocyclic double-bonded oxygen.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # We expect the carbon that bears the carbonyl to be a carbon.
            if atom.GetAtomicNum() != 6:
                continue

            # Fetch neighbors that are not in the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check if the bond to neighbor is a double bond.
                bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
                    continue
                
                # Check that the neighbor is an oxygen with no extra heavy connections.
                if nbr.GetAtomicNum() == 8 and nbr.GetDegree() == 1:
                    # Found a candidate: the ring is an oxolane with an exocyclic carbonyl.
                    # Further improvement could be added (connectivity checks, ensuring the oxygen is not shared by another functionality),
                    # but this already filters many false positives.
                    return True, ("Found a 5-membered ring (oxolane) with one oxygen and four carbons, "
                                  "and a ring carbon bearing an exocyclic double-bonded oxygen.")
                    
    return False, "No tetrahydrofuranone substructure found"
    
# Example usage (if run as a script):
if __name__ == "__main__":
    test_smiles = "CCCC(=O)N[C@@H]1CCOC1=O"  # Example: N-butyryl-L-homoserine lactone
    result, reason = is_tetrahydrofuranone(test_smiles)
    print(result, reason)