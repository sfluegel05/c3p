"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: Glycerophosphoinositol
Definition: Any glycerophospholipid having the polar alcohol inositol esterified 
to the phosphate group at the sn-3 position of the glycerol backbone.
This program uses heuristics by looking for a glycerol backbone, an inositol head group,
and a phosphate atom that connects them.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    
    The algorithm checks that the molecule:
      1. Has a glycerol backbone (HOCH2-CHOH-CH2OH pattern).
      2. Contains an inositol ring (approximated by myo-inositol; a cyclohexane
         with hydroxyl groups on all ring carbons).
      3. Contains at least one phosphorus whose oxygen neighbors connect (via single bonds)
         to atoms in both the glycerol fragment and the inositol fragment.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as glycerophosphoinositol, False otherwise.
        str: Reason for the classification decision.
    """
    # Convert the SMILES to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a query for a glycerol backbone fragment.
    # Note: HOCH2-CHOH-CH2OH can be represented as OCC(O)CO.
    glycerol_query = Chem.MolFromSmarts("OCC(O)CO")
    if glycerol_query is None:
        return False, "Error creating glycerol query"
    
    if not mol.HasSubstructMatch(glycerol_query):
        return False, "No glycerol backbone found"

    # Define a query for myo-inositol.
    # The structure of myo-inositol is approximated by a cyclohexane with hydroxyls on every ring carbon.
    # Here we use a non-stereospecific query: C1(C(C(C(C(C1O)O)O)O)O)O
    inositol_query = Chem.MolFromSmiles("C1(C(C(C(C(C1O)O)O)O)O)O")
    if inositol_query is None:
        return False, "Error creating inositol query"

    if not mol.HasSubstructMatch(inositol_query):
        return False, "No inositol head-group found"

    # Get all matches (atom indices) for the glycerol and inositol fragments in the molecule.
    glycerol_matches = mol.GetSubstructMatches(glycerol_query)
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    
    # Create sets of atom indices that belong to any glycerol or inositol fragment match.
    glycerol_atoms = set()
    for match in glycerol_matches:
        glycerol_atoms.update(match)
    
    inositol_atoms = set()
    for match in inositol_matches:
        inositol_atoms.update(match)

    # Now, for glycerophosphoinositol we expect a phosphate group linking these two.
    # We look for a phosphorus atom (atomic num 15) that is connected to oxygen atoms.
    phosphorus_found = False
    connection_found = False  # phosphate that connects inositol and glycerol

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue  # not phosphorus
        phosphorus_found = True
        neighbors = atom.GetNeighbors()
        # We expect at least two oxygen neighbors (one going to inositol and one to glycerol)
        oxygen_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8]
        if len(oxygen_neighbors) < 2:
            continue  # too few connections; try next phosphorus

        connected_to_inositol = False
        connected_to_glycerol = False

        # For each oxygen neighbor, check if it is attached to any atom in inositol or glycerol set.
        for oxy in oxygen_neighbors:
            for nbr in oxy.GetNeighbors():
                # Skip the phosphorus atom
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetIdx() in inositol_atoms:
                    connected_to_inositol = True
                if nbr.GetIdx() in glycerol_atoms:
                    connected_to_glycerol = True

        if connected_to_inositol and connected_to_glycerol:
            connection_found = True
            break  # found a suitable phosphate

    if not phosphorus_found:
        return False, "No phosphorus (phosphate group) found"
    if not connection_found:
        return False, "No phosphate found that connects both the glycerol backbone and the inositol head-group"
    
    return True, "Molecule contains a glycerol backbone with a phosphate group that connects to an inositol head-group"

# Example usage:
if __name__ == "__main__":
    # You can test the function with one of the example SMILES strings
    test_smiles = "P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@@H](O)COC(=O)CCCCCCCCCCCCCC)(O)=O"
    result, reason = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)