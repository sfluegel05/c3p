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
    Determines if a molecule is a flavonol based on its SMILES string.
    Flavonols are defined as 2-phenylchromen-4-ones that have a hydroxy 
    substituent (or an O–glycoside) at position 3 (i.e. the ring hydrogen 
    of a flavone replaced by –OH).
    
    The function first checks for the flavone core with a mapped atom (atom map 3),
    which must carry an exocyclic oxygen (atom mapped as 4). After a substructure
    match is found, the oxygen substituent is examined:
      • If the O has no other neighbor (besides the ring carbon) then it is considered a free –OH.
      • If the only extra neighbor is a carbon that is CH3 (a methoxy group),
        then the flavone is not a flavonol.
      • Otherwise (e.g. glycosides or larger groups) the 3-O functionality is acceptable.
    
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

    # Define a SMARTS for a 2-phenylchromen-4-one core with mapping.
    # The query requires:
    #   • A benzene ring (the 2-phenyl substituent) attached via a single bond (the dash)
    #   • To a heterocycle: c2oc([c:3]([O:4])c2=O)
    #     where the atom mapped as "3" must carry an oxygen substituent (mapped "4")
    query_smarts = "c1ccc(-c2oc([c:3]([O:4])c2=O))cc1"
    query = Chem.MolFromSmarts(query_smarts)
    if query is None:
        return False, "Error in generating the SMARTS query"

    # Look for substructure matches.
    # Use GetSubstructMatches; note that multiple matching fragments might be found.
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Molecule does not contain a 3-hydroxyflavone core required for flavonols"

    # To know which atom in the molecule corresponds to the mapped atoms, build a mapping:
    # (In our SMARTS, atom with map num 3 is the heterocyclic carbon expected to be at C3,
    # and it is bonded to an oxygen with map num 4.)
    mapnum_to_queryidx = {}
    for i, atom in enumerate(query.GetAtoms()):
        mapnum = atom.GetAtomMapNum()
        if mapnum:
            mapnum_to_queryidx[mapnum] = i

    # Now, for each match, retrieve the atom indices for mappings 3 and 4.
    # Then check the oxygen substituent (atom mapped 4) to ensure it is not a methoxy.
    valid_match_found = False
    for match in matches:
        # match is a tuple of indices corresponding to query atoms in order.
        # Using our mapping dictionary, get the atom index for map 3 and map 4.
        try:
            idx3 = match[ mapnum_to_queryidx[3] ]
            idx4 = match[ mapnum_to_queryidx[4] ]
        except KeyError:
            continue  # skip if mapping not found (should not happen)

        atom3 = mol.GetAtomWithIdx(idx3)
        atom4 = mol.GetAtomWithIdx(idx4)
        # Confirm that atom4 is exocyclic: besides its bond to atom3 (the ring carbon),
        # examine its other neighbors.
        o_neighbors = [nbr for nbr in atom4.GetNeighbors() if nbr.GetIdx() != idx3]
        
        # If there are no extra neighbors, then the substituent is a free –OH.
        if not o_neighbors:
            valid_match_found = True
            break
        # If there is exactly one extra neighbor, check if this neighbor is a methyl group.
        if len(o_neighbors) == 1:
            nbr = o_neighbors[0]
            # For a methoxy group, the oxygen will be attached to a single carbon that
            # carries three (implicit) hydrogens. (Use GetTotalNumHs to count H's)
            if nbr.GetAtomicNum() == 6 and nbr.GetTotalNumHs() == 3:
                # This match indicates a methoxy at position 3,
                # which disqualifies the structure as a flavonol.
                continue
            else:
                valid_match_found = True
                break
        else:
            # If the oxygen bears a substituent other than just a methyl group (e.g. a sugar),
            # we will accept this match.
            valid_match_found = True
            break

    if valid_match_found:
        return True, "Molecule contains a valid 3-hydroxyflavone core (flavonol scaffold)"
    else:
        return False, "The oxygen substituent at position 3 appears to be a methoxy group rather than a hydroxy or glycoside function."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test example: tambulin (should be classified as a flavonol)
    test_smiles = "COc1ccc(cc1)-c1oc2c(OC)c(OC)cc(O)c2c(=O)c1O"
    classification, explanation = is_flavonols(test_smiles)
    print("SMILES:", test_smiles)
    print("Classification:", classification)
    print("Reason:", explanation)