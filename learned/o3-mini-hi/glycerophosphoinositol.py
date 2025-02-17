"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: Glycerophosphoinositol
Definition: Any glycerophospholipid having the polar alcohol inositol esterified 
to the phosphate group at the sn-3 position of the glycerol backbone.

Improvement: This version not only looks for a glycerol backbone (approximated as OCC(O)CO)
and an inositol moiety (approximated by a cyclohexane ring with hydroxyls) but also inspects
each phosphorus atom. We “clean” the oxygen neighbors so that we ignore phosphate groups that 
bridge additional phosphates (i.e. in multiphosphates) and require that one oxygen (via a single bond)
directly connects to an atom from the glycerol fragment and one oxygen connects to an atom 
from the inositol fragment.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    
    The algorithm checks that the molecule:
      1. Contains a glycerol backbone fragment. (We approximate this with the pattern OCC(O)CO.)
      2. Contains an inositol head-group (approximated by a cyclohexane decorated with hydroxyl groups).
      3. Contains at least one phosphorus atom (the phosphate) that (via “clean” oxygen bonds) is
         connected on one side to the glycerol fragment and on the other side to the inositol fragment.
         In our check we filter out oxygen atoms that are linked to more than one phosphorus – as 
         such situations often occur in di- or polyphosphorylated compounds which should not be classified.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as glycerophosphoinositol, False otherwise.
        str: Reason for the classification decision.
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Query to find a glycerol backbone fragment.
    # Here we use a simplified pattern for glycerol: HOCH2-CHOH-CH2OH (as OCC(O)CO).
    glycerol_query = Chem.MolFromSmarts("OCC(O)CO")
    if glycerol_query is None:
        return False, "Error creating glycerol query"
    if not mol.HasSubstructMatch(glycerol_query):
        return False, "No glycerol backbone found"
    
    # Query for inositol.
    # We approximate myo-inositol by a cyclohexane ring with hydroxyl groups on every carbon.
    inositol_query = Chem.MolFromSmiles("C1(C(C(C(C(C1O)O)O)O)O)O")
    if inositol_query is None:
        return False, "Error creating inositol query"
    if not mol.HasSubstructMatch(inositol_query):
        return False, "No inositol head-group found"
    
    # Get all matching atoms from glycerol and inositol fragments.
    glycerol_matches = mol.GetSubstructMatches(glycerol_query)
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    glycerol_atoms = set()
    for match in glycerol_matches:
        glycerol_atoms.update(match)
    inositol_atoms = set()
    for match in inositol_matches:
        inositol_atoms.update(match)
    
    # Check for a phosphorus atom that bridges glycerol and inositol.
    phosphate_link_found = False
    phosphorus_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        phosphorus_found = True
        neighbors = atom.GetNeighbors()
        # Get all oxygen neighbors.
        oxygen_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8]
        if len(oxygen_neighbors) < 2:
            continue  # Not enough oxygens attached to a typical phosphate
        
        # Filter oxygens that are not "clean" – i.e. ignore any oxygen that is also attached to another phosphorus.
        clean_oxygens = []
        for oxy in oxygen_neighbors:
            # Count phosphorus neighbors other than current atom
            other_p_count = sum(1 for nbr in oxy.GetNeighbors() if nbr.GetAtomicNum() == 15 and nbr.GetIdx() != atom.GetIdx())
            if other_p_count == 0:
                clean_oxygens.append(oxy)
        if len(clean_oxygens) < 2:
            continue  # Likely part of a polyphosphate group; skip
        
        # Now check if at least one clean oxygen connects (via its other neighbor) to glycerol 
        # and at least one (could be the same or different) connects to inositol.
        connects_glycerol = False
        connects_inositol = False
        for oxy in clean_oxygens:
            for nbr in oxy.GetNeighbors():
                # Skip the phosphorus itself.
                if nbr.GetIdx() == atom.GetIdx():
                    continue
                if nbr.GetIdx() in glycerol_atoms:
                    connects_glycerol = True
                if nbr.GetIdx() in inositol_atoms:
                    connects_inositol = True
            # If both connections are found, no need to check further.
            if connects_glycerol and connects_inositol:
                break
        
        if connects_glycerol and connects_inositol:
            phosphate_link_found = True
            break  # Found a candidate phosphate bridging glycerol and inositol
    
    if not phosphorus_found:
        return False, "No phosphorus (phosphate group) found"
    if not phosphate_link_found:
        return False, "No phosphate found that cleanly connects the glycerol backbone and inositol head-group"
    
    return True, "Molecule contains a glycerol backbone with a phosphate group that connects directly to an inositol head-group"

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = "P(OC1C(O)C(O)C(O)[C@@H](O)C1O)(OC[C@@H](O)COC(=O)CCCCCCCCCCCCCC)(O)=O"
    result, reason = is_glycerophosphoinositol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)