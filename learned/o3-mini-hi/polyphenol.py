"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: polyphenol
Defined as: members of the class of phenols that contain 2 or more benzene rings,
each of which is substituted by at least one hydroxy (–OH) group (or derivative such as an O‐glycoside).
We use a heuristic: for each six‐membered aromatic ring (benzene), we check whether
there exists at least one exocyclic oxygen attached via a single bond that is not part of an unsaturated (e.g. carbonyl) system.
This relaxed criterion will allow recognition of both free phenolic –OH groups as well as cases where the –OH is masked (e.g. as an O‐glycoside).
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on the following heuristic:
    • The molecule must have at least 2 benzene rings (6-membered aromatic rings whose atoms are all carbons).
    • For each such ring, at least one ring carbon must be substituted by an oxygen atom via a single bond,
      where that oxygen is not a carbonyl oxygen (i.e. not double‐bonded to another carbon).
      This relaxed rule accepts both free –OH groups (oxygen bearing an explicit hydrogen)
      and cases where the hydroxyl has been substituted with a glycoside.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a polyphenol, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES to get a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that substitutions are clearer.
    mol = Chem.AddHs(mol)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {e}"

    # Get ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    valid_benzene_ring_count = 0

    # Iterate over each detected ring.
    for ring in atom_rings:
        # Only consider six-membered rings.
        if len(ring) != 6:
            continue
        
        # Check that every atom in the ring is a carbon and is aromatic.
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue

        # For each candidate benzene ring, look for an oxygen substituent that meets our relaxed criteria.
        found_valid_substituent = False
        for atom_idx in ring:
            ring_atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in ring_atom.GetNeighbors():
                # Only consider neighbors not in the ring (exocyclic substituents).
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # Get the bond between the ring carbon and the oxygen.
                    bond = mol.GetBondBetweenAtoms(ring_atom.GetIdx(), nbr.GetIdx())
                    if bond is None:
                        continue
                    # Consider only single bonds.
                    if bond.GetBondTypeAsDouble() != 1.0:
                        continue

                    # Now check that this oxygen is not part of a carbonyl arrangement.
                    # In a carbonyl oxygen the oxygen is double-bonded to another carbon.
                    skip = False
                    for oxy_bond in nbr.GetBonds():
                        # Skip the bond to the ring atom.
                        if oxy_bond.GetBeginAtomIdx() == ring_atom.GetIdx() or oxy_bond.GetEndAtomIdx() == ring_atom.GetIdx():
                            continue
                        # If any other bond from oxygen is a double bond, we assume the oxygen is in a carbonyl.
                        if oxy_bond.GetBondTypeAsDouble() > 1.0:
                            skip = True
                            break
                    if skip:
                        continue

                    # Passed our criteria – this oxygen substituent is acceptable.
                    found_valid_substituent = True
                    break
            if found_valid_substituent:
                break
        
        if found_valid_substituent:
            valid_benzene_ring_count += 1

    if valid_benzene_ring_count >= 2:
        return True, (f"Found {valid_benzene_ring_count} benzene ring(s) each with a valid oxygen substituent "
                      "(not involved in a carbonyl), classifying as a polyphenol")
    else:
        return False, (f"Only {valid_benzene_ring_count} benzene ring(s) with a valid oxygen substituent were found")

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "O1C2=C(OC)C(=C(C3=C(O)C=4OCOC4C(=C3C)OC)C(=C2OC1)O)C",  # Benzocamphorin E
        "OC1C(O)c2c(O)cc(O)cc2OC1c1cc(O)c(O)c(O)c1",  # flavan-3,3',4,4',5,5',7-heptol
        "O=C(OC1=CC(O)=CC(=C1)CCCCC)C2=C(O[C@@H]3O[C@H](C(=O)O)[C@H](O)[C@@H]([C@H]3O)O)C=C(O)C=C2CCCCC",  # Ascotricin B
        # False negatives (previously missed; may now be captured due to relaxed criteria):
        "O([C@H]1[C@H](O)C(O[C@@H](OC2=CC=3OC=C(C(=O)C3C=C2)C4=CC=5OCOC5C=C4)C1O)CO)[C@@H]6OC([C@@H](O)[C@H](O)C6O)CO",  # Pseudobaptigenin 7-O-laminaribioside
    ]
    for smi in test_smiles:
        result, reason = is_polyphenol(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)