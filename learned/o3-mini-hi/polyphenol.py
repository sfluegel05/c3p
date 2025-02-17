"""
Classifies: CHEBI:26195 polyphenol
"""
"""
Classifies: polyphenol
Defined as 'Members of the class of phenols that contain 2 or more benzene rings
each of which is substituted by at least one free (non-glycosylated/non-sulfated) hydroxy group.'
In this implementation we attempt to count only those aromatic rings (6-membered, all aromatic C)
that have at least one directly attached hydroxyl group (–OH) that is not covalently linked to another heavy atom.
"""
from rdkit import Chem

def is_polyphenol(smiles: str):
    """
    Determines if a molecule is a polyphenol based on its SMILES string.
    We require the molecule to contain at least 2 aromatic rings (benzene rings) 
    and on each such ring at least one substituent which is a free hydroxy group (-OH).
    
    We define a benzene ring as a 6-membered ring where every atom is an aromatic carbon.
    For a hydroxy substituent, we require an oxygen atom connected to the ring carbon that
    has at least one attached hydrogen and is not also connected to another heavy atom (i.e. not part of a glycoside or sulfate).
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a polyphenol, False otherwise.
        str: Reason for the classification decision.
    """
    # parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # update cached properties (including aromaticity and ring info)
    mol.UpdatePropertyCache()
    
    # get ring info: each ring is represented as a tuple of atom indices
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    # count benzene rings with at least one “free” hydroxy group attached.
    hydroxy_ring_count = 0

    # iterate over all rings detected in the molecule
    for ring in atom_rings:
        # we are only looking for 6-membered rings
        if len(ring) != 6:
            continue

        # check if all atoms are aromatic carbons (atomic number 6 and aromatic)
        is_benzene = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue

        # look for a free hydroxy group (-OH) attached to any of the ring atoms
        has_free_OH = False
        for idx in ring:
            ring_atom = mol.GetAtomWithIdx(idx)
            for nbr in ring_atom.GetNeighbors():
                # if the neighbor is part of the ring, skip it
                if nbr.GetIdx() in ring:
                    continue
                # we are interested in an oxygen neighbor
                if nbr.GetAtomicNum() == 8:
                    # To qualify as a free OH, the oxygen should have at least one H,
                    # and it should not be attached to another heavy atom (atomic number > 1)
                    # besides the ring carbon.
                    # In glycosidic or ester bonds the oxygen is attached to a second heavy atom.
                    heavy_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() > 1]
                    # We expect a free hydroxy oxygen to be attached to exactly 1 heavy atom (the ring carbon)
                    if len(heavy_neighbors) == 1:
                        # also check that it carries at least one hydrogen (explicit + implicit)
                        if nbr.GetTotalNumHs() >= 1:
                            has_free_OH = True
                            break
            if has_free_OH:
                break

        if has_free_OH:
            hydroxy_ring_count += 1

    # Decision based on the number of benzene rings with a free -OH group.
    if hydroxy_ring_count < 2:
        return False, f"Only {hydroxy_ring_count} benzene ring(s) with a free hydroxy group found"
    else:
        return True, f"Found {hydroxy_ring_count} benzene rings each substituted with a free hydroxy group, classifying as a polyphenol"

# For example testing purpose (you may uncomment and run locally):
# if __name__ == "__main__":
#     test_smiles = "OC1C(O)c2c(O)cc(O)cc2OC1c1cc(O)c(O)c(O)c1"  # example: flavan-3,3',4,4',5,5',7-heptol
#     result, reason = is_polyphenol(test_smiles)
#     print("Result:", result)
#     print("Reason:", reason)