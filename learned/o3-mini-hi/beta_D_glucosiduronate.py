"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate
Definition:
  A carbohydrate acid derivative anion obtained by deprotonation of the carboxy group
  of any beta-D-glucuronic acid; major species at pH 7.3.
  
This function checks for a glucuronate-like fragment by searching for a pyranose (6-membered) ring
with exactly one ring oxygen and the presence of an exocyclic carboxylate group (C(=O)[O-])
attached to one of the ring carbons.
"""

from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule has a beta-D-glucosiduronate fragment.
    This simplified algorithm looks for:
      - A valid molecule.
      - At least one six-membered ring (pyranose-like) containing five carbons and one oxygen.
      - One of the ring carbons (attained from the CH2 portion) attached exocyclically to a carboxylate group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule contains a beta-D-glucosiduronate fragment, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for a carboxylate group.
    carboxylate_smarts = "C(=O)[O-]"
    carboxylate_query = Chem.MolFromSmarts(carboxylate_smarts)
    if carboxylate_query is None:
        # Should not happen.
        return False, "Carboxylate SMARTS pattern error"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings in molecule"
    
    # Helper: Check if an atom (carboxylate candidate) is indeed a carboxylate carbon.
    def is_carboxylate_carbon(atom):
        # Must be a carbon.
        if atom.GetAtomicNum() != 6:
            return False
        dbl_oxygen = False
        single_oxygen = False
        # Loop over neighbors.
        for nbr in atom.GetNeighbors():
            # Get bond info.
            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
            if nbr.GetAtomicNum() == 8:
                # Check bond order.
                if bond.GetBondTypeAsDouble() == 2:
                    dbl_oxygen = True
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    # Check if oxygen bears a negative formal charge.
                    if nbr.GetFormalCharge() == -1:
                        single_oxygen = True
        return dbl_oxygen and single_oxygen

    # Loop through each ring.
    for ring_atoms in ring_info.AtomRings():
        if len(ring_atoms) != 6:
            continue  # we only consider six-membered rings (pyranose rings)
        # Count ring oxygen atoms.
        oxy_in_ring = 0
        carbon_indices = []
        for idx in ring_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                oxy_in_ring += 1
            elif atom.GetAtomicNum() == 6:
                carbon_indices.append(idx)
        # A pyranose ring typically has one oxygen and five carbons.
        if oxy_in_ring != 1 or len(carbon_indices) != 5:
            continue
        
        # Now, check if any ring carbon has an exocyclic neighbor that is a carboxylate group.
        for c_idx in carbon_indices:
            ring_atom = mol.GetAtomWithIdx(c_idx)
            # Look at neighbors that are NOT in the ring.
            for nbr in ring_atom.GetNeighbors():
                if nbr.GetIdx() in ring_atoms:
                    continue
                # If the neighbor is carbon and may be a carboxylate carbon.
                if nbr.GetAtomicNum() == 6 and is_carboxylate_carbon(nbr):
                    # Additionally, one can check that the overall molecule has the carboxylate SMARTS.
                    if mol.HasSubstructMatch(carboxylate_query):
                        # Optionally, one might want to check for additional hydroxy groups on the ring carbons,
                        # but here we assume the presence of the carboxylate on a pyranose ring is sufficient.
                        return True, "Found a six-membered pyranose ring with an exocyclic carboxylate group consistent with beta-D-glucosiduronate"
        
    return False, "No beta-D-glucosiduronate fragment found"
    
# Example usage:
if __name__ == "__main__":
    example_smiles = "O[C@@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]"
    result, reason = is_beta_D_glucosiduronate(example_smiles)
    print("Result:", result)
    print("Reason:", reason)