"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: 2-hydroxy fatty acid

Definition:
A 2-hydroxy fatty acid is an acyclic (non‐ring) aliphatic molecule that bears
exactly one terminal carboxylic acid group. In addition, the unique carbon attached
to that acid group (the “α–carbon”) must carry exactly one hydroxyl (-OH) substituent.
No other oxygens (whether –OH or carbonyl) are allowed in the molecule. (Any extra
stereochemistry in the alkyl chain, if present, is ignored for classification.)
"""

from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    Criteria:
      1) The molecule must be acyclic (no rings).
      2) It must have exactly one terminal carboxylic acid group.
         (Identified by the pattern: "C(=O)[O;H1]". The acid carbon is expected to be
          bonded to two oxygens (one as a carbonyl and one as hydroxyl) and exactly one carbon.)
      3) The unique (α–) carbon attached to the acid carbon must have exactly one hydroxyl (-OH)
         substituent (its only oxygen neighbor other than the acid carbon).
      4) Aside from the acid group (acid carbon with its two oxygens) and the alpha –OH, no additional
         oxygen atoms (or oxygen functionality) may be present.
      5) The non–acid part (the fatty chain) must be composed almost entirely of carbon atoms.
      6) Extra explicit chiral centers (other than possibly the alpha–carbon) will be ignored so as not
         to reject branched but valid fatty acids.
      7) The fatty chain length (all carbons except the acid carbon) should be at least 3.
         
    Args:
         smiles (str): SMILES string of the molecule.
         
    Returns:
         bool: True if the molecule qualifies as a 2-hydroxy fatty acid, else False.
         str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Reject molecules with rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), not a typical fatty acid"

    # 2. Identify the terminal carboxylic acid group.
    # Use SMARTS to find a carboxyl group.
    # Note: The SMARTS "C(=O)[O;H1]" looks for a carbon doubly bonded to O and singly bonded to an -OH.
    carboxyl_smarts = "C(=O)[O;H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    # We expect exactly one acid carbon. In each match, atom 0 of the SMARTS is the acid carbon.
    acid_carbon_idxs = set(match[0] for match in carboxyl_matches)
    if len(acid_carbon_idxs) != 1:
        return False, "Molecule must have exactly one terminal carboxylic acid group"
    acid_carbon_idx = list(acid_carbon_idxs)[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)

    # Verify that the acid carbon is terminal:
    # It should be bonded to exactly: one carbon neighbor, one oxygen involved in the carbonyl,
    # and one oxygen in the -OH part.
    acid_neighbors = acid_carbon.GetNeighbors()
    oxy_neighbors = [nbr for nbr in acid_neighbors if nbr.GetAtomicNum() == 8]
    carbon_neighbors = [nbr for nbr in acid_neighbors if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1 or len(oxy_neighbors) != 2:
        return False, "Acid carbon is not terminal (unexpected number of neighbors)"

    # Record the indices of oxygens of the carboxyl group.
    allowed_oxygens = {nbr.GetIdx() for nbr in oxy_neighbors}

    # 3. Identify the α–carbon: the unique carbon neighbor of the acid carbon.
    alpha_carbon = carbon_neighbors[0]

    # 4. Check that the α–carbon has exactly one hydroxyl (-OH) substituent.
    # Look among the neighbors of the α–carbon (ignoring the acid carbon)
    alpha_oh_atoms = []
    for nbr in alpha_carbon.GetNeighbors():
        if nbr.GetIdx() == acid_carbon_idx:
            continue
        if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
            alpha_oh_atoms.append(nbr)
    if len(alpha_oh_atoms) != 1:
        return False, f"Alpha carbon must have exactly one hydroxyl substituent, found {len(alpha_oh_atoms)}"
    alpha_oh = alpha_oh_atoms[0]
    allowed_oxygens.add(alpha_oh.GetIdx())  # add the alpha –OH oxygen to allowed set

    # 5. Now ensure that no other oxygen is present in the molecule.
    # (Any oxygen not in the allowed set would represent an extra functional group.)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            if atom.GetIdx() not in allowed_oxygens:
                return False, f"Extra oxygen functionality found at atom index {atom.GetIdx()}"

    # 6. Check that the fatty (non–acid) chain is aliphatic.
    # Exclude the acid carbon and its two oxygen neighbors.
    excluded_idxs = set([acid_carbon_idx])
    for idx in allowed_oxygens:
        excluded_idxs.add(idx)
    non_acid_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() not in excluded_idxs]
    if not non_acid_atoms:
        return False, "No carbon chain found outside the acid group"
    n_nonacid = len(non_acid_atoms)
    n_carbons = sum(1 for atom in non_acid_atoms if atom.GetAtomicNum() == 6)
    # We require that nearly all of the rest are carbons.
    if n_carbons / n_nonacid < 0.75:
        return False, f"Non-carboxyl part is not aliphatic enough (carbon fraction = {n_carbons / n_nonacid:.2f})"

    # 7. Check the chain length.
    # Count total carbons and subtract the acid carbon.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    chain_length = total_carbons - 1  # subtract the acid carbon
    if chain_length < 3:
        return False, f"Carbon chain too short ({chain_length} carbons) for a fatty acid"

    # 8. (Chirality considerations)
    # Although the definition mentions that if chiral centers are defined the only allowed should be the alpha-carbon,
    # many natural fatty acids have extra chiral atoms with undefined stereochemistry.
    # Thus, we choose to ignore chiral centers beyond the alpha carbon.
    # (If needed one could refine this check by allowing only undefined chiral centers outside the alpha carbon.)

    reason = (f"Found a terminal carboxylic acid group (acid carbon idx {acid_carbon_idx}) with "
              f"an alpha carbon (idx {alpha_carbon.GetIdx()}) bearing a unique hydroxyl substituent, "
              f"and a predominantly aliphatic chain (chain length = {chain_length} carbons, "
              f"carbon fraction = {n_carbons/n_nonacid:.2f}), classifying the molecule as a 2-hydroxy fatty acid.")
    return True, reason


# Example usage:
if __name__ == "__main__":
    # Test one of the true positive examples:
    test_smiles = "CCCCCCCCCCCC[C@@H](O)C(O)=O"  # (2R)-2-hydroxytetradecanoic acid
    result, reason = is_2_hydroxy_fatty_acid(test_smiles)
    print(result, reason)