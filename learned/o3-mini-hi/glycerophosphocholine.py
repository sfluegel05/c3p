"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
Definition: The glycerol phosphate ester of a phosphocholine.
A true glycerophosphocholine (GPC) should contain (i) a phosphocholine headgroup – an oxygen‐bonded 
ethylene trimethylammonium fragment attached to a phosphate group – and (ii) that phosphate must be
linked via an oxygen to a primary carbon (typically from a glycerol or glycerol‐like moiety). 
This heuristic is meant to help distinguish true glycerophosphocholines from related phospholipids 
(such as sphingomyelins or diacyl phosphatidylcholines) that contain a phosphocholine fragment but have 
a different backbone.
Note: This is only an heuristic method.
"""

from rdkit import Chem

def is_glycerophosphocholine(smiles: str):
    """
    Determines whether a molecule (given as a SMILES string) is a glycerophosphocholine.
    This implementation proceeds roughly in two steps:
    
      1. It verifies that there is at least one phosphorus atom and at least one phosphocholine‐like 
         substructure. We search for a substructure match that contains the fragment:
           CH2–O–P(=O)(O or [O-])–OCC[N+](C)(C)C
         (allowing small variations in the protonation state).
    
      2. For at least one phosphorus atom having such a headgroup, we inspect one of its oxygen neighbors 
         (the one connecting to the headgroup branch). We require that the oxygen connects to a carbon that 
         is “primary” (which here means it has (at least) two implicit hydrogens). This is our rough proxy 
         for a glycerol (or glycerol‐like) connectivity.
         
    Args:
      smiles (str): SMILES string for the molecule.
      
    Returns:
      bool: True if the molecule is accepted as glycerophosphocholine, else False.
      str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1. Check that at least one phosphorus atom is present.
    if not any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "No phosphorus atom found; not a phospholipid"
    
    # Define a list of SMARTS patterns for a phosphocholine headgroup that also enforce a primary CH2 at 
    # the connection to the phosphate.
    # We require that an oxygen attached to the phosphorus is in turn attached to a CH2 group.
    headgroup_smarts_list = [
        "[CH2]O[P](=O)(O)OCC[N+](C)(C)C",      # all oxygens protonated
        "[CH2]O[P](=O)([O-])OCC[N+](C)(C)C",    # one deprotonated oxygen (pattern ordering variant 1)
        "[CH2]O[P]([O-])(=O)OCC[N+](C)(C)C"     # one deprotonated oxygen (pattern ordering variant 2)
    ]
    
    # Look for at least one match of any of these patterns.
    headgroup_match_found = False
    for smarts in headgroup_smarts_list:
        patt = Chem.MolFromSmarts(smarts)
        if patt is None:
            continue
        if mol.HasSubstructMatch(patt):
            headgroup_match_found = True
            break
    if not headgroup_match_found:
        return False, "Phosphocholine headgroup substructure not found"
    
    # Step 2. Look for evidence of proper (glycerol-like) connectivity.
    # For each phosphorus atom, look at its oxygen neighbors.
    # We expect that one oxygen (not the one carrying the choline branch) is attached to a primary carbon.
    glycerol_connectivity_found = False
    for p_atom in mol.GetAtoms():
        if p_atom.GetAtomicNum() != 15:
            continue
        for o_atom in p_atom.GetNeighbors():
            # Consider only oxygen neighbors.
            if o_atom.GetAtomicNum() != 8:
                continue
            # Get the neighbor(s) of this oxygen besides the phosphorus.
            other_neighbors = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetIdx() != p_atom.GetIdx()]
            # Skip if an oxygen neighbor is directly part of the choline branch (look for N+ among neighbors).
            choline_branch = False
            for nbr in other_neighbors:
                if nbr.GetAtomicNum() == 7 and nbr.GetFormalCharge() == 1:
                    choline_branch = True
                    break
            if choline_branch:
                continue
            # Now, for each remaining neighbor that is a carbon, check if it appears to be primary.
            for nbr in other_neighbors:
                if nbr.GetAtomicNum() == 6:
                    # Count the (implicit+explicit) hydrogens on the carbon.
                    nHs = nbr.GetTotalNumHs()
                    # For a CH2 group, typically expect two hydrogens.
                    if nHs >= 2:
                        glycerol_connectivity_found = True
                        break
            if glycerol_connectivity_found:
                break
        if glycerol_connectivity_found:
            break
    if not glycerol_connectivity_found:
        return False, "Glycerol backbone connectivity not found (phosphate not linked via an oxygen to a primary CH2 carbon)"
    
    return True, "Molecule contains a phosphocholine headgroup with proper glycerol connectivity"

# Example usage.
if __name__ == "__main__":
    # Test with one of the provided SMILES strings (example: 2-[(9Z)-12-hydroxyoctadec-9-enoyl]-sn-glycero-3-phosphocholine).
    test_smiles = "P(OC[C@@H](CO)OC(CCCCCCC/C=C\\CC(CCCCCC)O)=O)(=O)(OCC[N+](C)(C)C)[O-]"
    result, reason = is_glycerophosphocholine(test_smiles)
    print(result, reason)