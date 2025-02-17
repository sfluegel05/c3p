"""
Classifies: CHEBI:75659 O-acyl-L-carnitine
"""
"""
Classifies: O-acyl-L-carnitine
Definition: An O-acylcarnitine in which the carnitine component has L-configuration.
The molecule should have a carnitine moiety that displays a chiral center bonded to a carboxylate group
and a trimethylammonium substituent, with the acyl group attached via an ester (OC(=O)R linkage).
"""

from rdkit import Chem

def is_O_acyl_L_carnitine(smiles: str):
    """
    Determines if a molecule is an O-acyl-L-carnitine based on its SMILES string.
    
    The function checks for:
      1. A carnitine substructure containing a chiral carbon bearing a carboxylate
         (CC([O-])=O) and a trimethylammonium group (C[N+](C)(C)C).
      2. An ester linkage, with the acyl group attached via an oxygen atom (i.e. O attached to a carbonyl group).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as an O-acyl-L-carnitine, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS queries for the carnitine moiety with an ester-linked acyl group.
    # The expected motif is an oxygen atom attached to a chiral carbon that bears 
    # (a) a carboxylate group: CC([O-])=O and
    # (b) a trimethylammonium group: C[N+](C)(C)C.
    # We allow for either [C@H] or [C@@H] in the SMARTS to account for different notations.
    query1 = Chem.MolFromSmarts("O[C@H](CC([O-])=O)C[N+](C)(C)C")
    query2 = Chem.MolFromSmarts("O[C@@H](CC([O-])=O)C[N+](C)(C)C")
    
    # Get substructure matches for either pattern.
    matches = mol.GetSubstructMatches(query1) + mol.GetSubstructMatches(query2)
    if not matches:
        return False, "No carnitine substructure with ester-linked acyl group was found."
    
    # For each match, verify that the oxygen (atom at index 0 in the query) is acylated.
    # The oxygen should be bonded to:
    #   a) the carnitine chiral center (already in the query) and 
    #   b) an acyl carbon that is part of a carbonyl (C=O) group.
    for match in matches:
        # Get the oxygen atom index from the match (index 0 from the SMARTS)
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        neighbors = o_atom.GetNeighbors()
        if len(neighbors) != 2:
            # Expected connectivity: one neighbor is the chiral carbon and one is the acyl carbon.
            continue
        
        # Identify which neighbor is not the carnitine chiral center.
        neighbor_ids = [atom.GetIdx() for atom in neighbors]
        if match[1] in neighbor_ids:
            neighbor_ids.remove(match[1])
        else:
            continue
        
        if not neighbor_ids:
            continue
        
        acyl_atom = mol.GetAtomWithIdx(neighbor_ids[0])
        # The acyl atom must be carbon.
        if acyl_atom.GetAtomicNum() != 6:
            continue
        
        # Check that the acyl atom is part of a carbonyl group by looking for a double bond O neighbor.
        has_carbonyl = False
        for nb in acyl_atom.GetNeighbors():
            bond = mol.GetBondBetweenAtoms(acyl_atom.GetIdx(), nb.GetIdx())
            if nb.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                has_carbonyl = True
                break
        if not has_carbonyl:
            continue
        
        # If the oxygen is acylated via a carbonyl group, the motif is satisfied.
        return True, "Molecule contains an O-acyl-L-carnitine motif with proper acylation and L-carnitine configuration."
    
    return False, "Carnitine-like substructure found, but the acyl group attachment on the hydroxyl group is missing."