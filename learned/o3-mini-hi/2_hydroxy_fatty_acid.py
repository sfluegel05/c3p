"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: 2-hydroxy fatty acid

Definition:
A 2-hydroxy fatty acid is an acyclic (non‐ring) aliphatic molecule 
that bears exactly one terminal carboxylic acid group. In addition, 
the unique carbon attached to that acid group (the “α–carbon”) must 
carry exactly one hydroxyl (-OH) substituent and no additional –OH groups 
should be present in the carbon chain. Finally, the “fatty” chain must be 
predominantly composed of carbon atoms and (if chiral) only the α–carbon is allowed 
to have defined stereochemistry.
"""

from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    Criteria:
      1) The molecule is acyclic.
      2) It has exactly one terminal carboxylic acid group.
         We look for the pattern "C(=O)[O;H1]" (the acid carbon).
         Additionally we require that the acid carbon be terminal
         (having exactly one carbon neighbor).
      3) The unique (α–) carbon attached to the acid carbon must have exactly one 
         hydroxyl (-OH) substituent.
      4) No other –OH groups (hydroxyls) are allowed on carbons in the molecule.
      5) The remainder of the molecule must be “aliphatic” –
         i.e. most non-acid atoms should be carbons.
      6) If any chiral centers are found, then (if present) the only chiral center 
         allowed is the α–carbon.
         
    Args:
         smiles (str): SMILES string of the molecule.
         
    Returns:
         bool: True if the molecule qualifies as a 2-hydroxy fatty acid, else False.
         str: A reason explaining the classification decision.
    """

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Reject if the molecule contains rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a typical fatty acid"
    
    # 2. Find the terminal carboxylic acid group.
    # SMARTS for carboxylic acid: carbon doubly bonded to O and singly bonded to an O with one H.
    carboxyl_smarts = "C(=O)[O;H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"
    
    # Get the acid carbon indices (the “C” in the SMARTS is atom 0 in each match).
    acid_carbons = set(match[0] for match in carboxyl_matches)
    if len(acid_carbons) != 1:
        return False, "Molecule must have exactly one terminal carboxylic acid group"
    
    acid_carbon_idx = list(acid_carbons)[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Verify that the acid carbon is terminal:
    # It should be bonded only to (i) one oxygen by a double bond, (ii) one hydroxyl oxygen, and (iii) one carbon.
    acid_neighbors = acid_carbon.GetNeighbors()
    oxygen_neighbors = [nbr for nbr in acid_neighbors if nbr.GetAtomicNum() == 8]
    carbon_neighbors = [nbr for nbr in acid_neighbors if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1 or len(oxygen_neighbors) != 2:
        return False, "Acid carbon is not terminal (unexpected number of neighbors)"
    
    # 3. Identify the α–carbon:
    alpha_carbon = carbon_neighbors[0]
    
    # 4. Check for hydroxyl substituent(s) on the α–carbon.
    # We look among the neighbors of the α–carbon excluding the acid carbon.
    alpha_hydroxyl = None
    for nbr in alpha_carbon.GetNeighbors():
        if nbr.GetIdx() == acid_carbon_idx:
            continue
        # We consider an oxygen atom with at least one hydrogen (a -OH).
        if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
            alpha_hydroxyl = nbr
            break
    if alpha_hydroxyl is None:
        return False, "Alpha carbon does not have a hydroxyl (-OH) substituent"
    
    # 5. Ensure no extra hydroxyl groups occur on any carbon in the molecule (i.e. the α–OH is unique).
    extra_oh_count = 0
    extra_oh_locations = []
    # For every oxygen in the molecule that is attached to a carbon and bears at least one hydrogen,
    # ignore those that are attached to the acid carbon (the carboxyl group).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 8:
            continue
        if atom.GetTotalNumHs() < 1:
            continue
        # Check if this oxygen is part of the carboxyl group:
        neighbor_ids = [nbr.GetIdx() for nbr in atom.GetNeighbors()]
        if acid_carbon_idx in neighbor_ids:
            continue  # skip carboxyl oxygens
        # Otherwise, count it as a hydroxyl group on the carbon chain.
        extra_oh_count += 1
        # Record the neighbor carbon(s)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6:
                extra_oh_locations.append(nbr.GetIdx())
    # We expect exactly one external -OH and it must be on the α–carbon.
    if extra_oh_count != 1 or alpha_carbon.GetIdx() not in extra_oh_locations:
        return False, f"Found {extra_oh_count} extra hydroxyl groups; expected exactly one on the alpha carbon"
    
    # 6. Check aliphatic character:
    # To be a fatty acid the remainder of the molecule (excluding the acid group)
    # should be composed almost entirely of carbons.
    # We will exclude the acid carbon and its two carboxyl oxygens.
    excluded_idxs = {acid_carbon_idx}
    for nbr in acid_carbon.GetNeighbors():
        if nbr.GetAtomicNum() == 8:
            excluded_idxs.add(nbr.GetIdx())
    non_acid_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() not in excluded_idxs]
    if non_acid_atoms:
        n_nonacid = len(non_acid_atoms)
        n_carbons = sum(1 for atom in non_acid_atoms if atom.GetAtomicNum() == 6)
        frac_c = n_carbons / n_nonacid
        # Set a threshold: most of these atoms should be carbons (at least 75%).
        if frac_c < 0.75:
            return False, f"Non-carboxyl part is not aliphatic enough (carbon fraction = {frac_c:.2f})"
    else:
        return False, "No non-carboxyl atoms found"
    
    # 7. Check stereochemistry: if there are any chiral centers, then (if present)
    # the only chiral center allowed should be the α–carbon.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    # chiral_centers is a list of tuples (atom_index, chirality).
    # If more than one chiral center is present, and it is not solely the α–carbon, then reject.
    if len(chiral_centers) > 1:
        # Allow cases with no chiral info or only the alpha carbon chiral
        if not all(center[0] == alpha_carbon.GetIdx() for center in chiral_centers):
            return False, f"Found extra chiral center(s): {chiral_centers}, not allowed for a typical 2-hydroxy fatty acid"
    
    # 8. (Optional) Check that there is a minimal aliphatic chain length.
    # Here we count all carbons not in the acid group.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    chain_length = total_carbons - 1  # subtract the acid carbon
    # We allow very short unsaturated acids (e.g. penta-2,4-dienoic acid) but for saturated systems
    # we enforce a slightly longer chain. (Many true examples have chain_length >= 4.)
    if chain_length < 3:
        return False, f"Carbon chain too short ({chain_length} carbons) for a fatty acid"
    
    return True, (f"Found a terminal carboxylic acid with an alpha carbon bearing a unique hydroxyl substituent, "
                  f"and an aliphatic chain (chain length = {chain_length} carbons, carbon fraction = {frac_c:.2f}), "
                  "classifying the molecule as a 2-hydroxy fatty acid")


# Example usage:
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCC[C@@H](O)C(O)=O"  # (2R)-2-hydroxytetradecanoic acid
    result, reason = is_2_hydroxy_fatty_acid(test_smiles)
    print(result, reason)