"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: 2-hydroxy fatty acid

Definition:
A 2–hydroxy fatty acid is an acyclic (non–ring) aliphatic molecule that bears
exactly one terminal carboxylic acid group. In addition, the unique carbon attached
to that acid group (the “alpha–carbon”) must carry exactly one hydroxyl (-OH)
substituent. No extra oxygen functionalities (beyond the two oxygens in the acid group
and the single alpha–OH) are allowed. Finally, a saturated fatty acid must have a sufficiently
long carbon chain (non–acid part ≥8 carbons) whereas unsaturated variants (showing C=C)
are allowed to be shorter.
"""

from rdkit import Chem

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as a 2–hydroxy fatty acid given its SMILES.
    
    The criteria are:
      1) The molecule must be acyclic (i.e. contain no rings).
      2) It must contain exactly one terminal carboxylic acid group (using the SMARTS
         "C(=O)[O;H1]").
      3) The unique (alpha–) carbon attached to the acid carbon must have exactly one
         hydroxyl (-OH) substituent (i.e. one oxygen attached that is hydrogenated).
      4) No other oxygen atoms should be present in the molecule.
      5) The molecule must contain only C, H, and O atoms.
      6) The fatty (non–acid) chain should be predominantly aliphatic.
      7) If the molecule is fully saturated (i.e. no C=C bonds) then the carbon chain
         (all carbons excluding the acid carbon) must be at least 8 atoms long.
         Unsaturated acids (with one or more C=C) are allowed to have shorter chains.
         
    Args:
         smiles (str): SMILES representation of the molecule.
         
    Returns:
         (bool, str): Tuple where the boolean indicates classification as a 2–hydroxy fatty acid,
                      and the reason details the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Reject molecules with rings.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring(s), not an acyclic fatty acid"

    # 2. Ensure that only allowed elements (C, H, O) are present.
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in {"C", "H", "O"}:
            return False, f"Atom {atom.GetSymbol()} is not allowed in a typical fatty acid"

    # 3. Identify the terminal carboxylic acid group.
    # Use SMARTS: acid group pattern "C(=O)[O;H1]".
    carboxyl_smarts = "C(=O)[O;H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if not carboxyl_matches:
        return False, "No carboxylic acid group found"

    # Expect exactly one unique acid carbon (atom 0 in each match).
    acid_carbon_idxs = set(match[0] for match in carboxyl_matches)
    if len(acid_carbon_idxs) != 1:
        return False, "Molecule must have exactly one terminal carboxylic acid group"
    acid_carbon_idx = list(acid_carbon_idxs)[0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)

    # Verify the acid carbon is terminal: one carbon neighbor and two oxygen neighbors.
    acid_neighbors = acid_carbon.GetNeighbors()
    oxy_neighbors = [nbr for nbr in acid_neighbors if nbr.GetAtomicNum() == 8]
    carbon_neighbors = [nbr for nbr in acid_neighbors if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1 or len(oxy_neighbors) != 2:
        return False, "Acid carbon is not terminal (unexpected bonding pattern)"

    # Allowed oxygens: those in the carboxyl group.
    allowed_oxygens = {nbr.GetIdx() for nbr in oxy_neighbors}

    # 4. Identify the alpha carbon (the sole carbon neighbor of the acid carbon).
    alpha_carbon = carbon_neighbors[0]

    # 5. Check alpha carbon has exactly one hydroxyl (-OH) substituent.
    alpha_oh_atoms = []
    for nbr in alpha_carbon.GetNeighbors():
        # Skip the acid carbon
        if nbr.GetIdx() == acid_carbon_idx:
            continue
        # Check if neighbor is oxygen and has at least one hydrogen
        if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
            alpha_oh_atoms.append(nbr)
    if len(alpha_oh_atoms) != 1:
        return False, f"Alpha carbon must have exactly one hydroxyl substituent; found {len(alpha_oh_atoms)}"
    alpha_oh = alpha_oh_atoms[0]
    allowed_oxygens.add(alpha_oh.GetIdx())

    # 6. Ensure no unexpected oxygen is present.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:
            if atom.GetIdx() not in allowed_oxygens:
                return False, f"Extra oxygen functionality found at atom index {atom.GetIdx()}"

    # 7. Check that the remainder of the molecule (non–acid part) is predominantly aliphatic.
    # Exclude the acid carbon and its two oxygens.
    excluded_idxs = {acid_carbon_idx} | allowed_oxygens
    non_acid_atoms = [atom for atom in mol.GetAtoms() if atom.GetIdx() not in excluded_idxs]
    if not non_acid_atoms:
        return False, "No carbon chain found outside the acid group"
    n_nonacid = len(non_acid_atoms)
    n_carbons = sum(1 for atom in non_acid_atoms if atom.GetAtomicNum() == 6)
    if n_carbons / n_nonacid < 0.75:
        return False, f"Non–acid part is not aliphatic enough (carbon fraction = {n_carbons/n_nonacid:.2f})"
    
    # 8. Count total number of carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Define the fatty chain length as all carbons except the acid carbon.
    chain_length = total_carbons - 1

    # 9. If molecule is saturated (no C=C), require a minimum chain length.
    # Note: unsaturated acids may be naturally shorter.
    unsaturation = mol.HasSubstructMatch(Chem.MolFromSmarts("C=C"))
    if not unsaturation:
        # For saturated fatty acids, require that the chain (non–acid part) has at least 8 carbons.
        if chain_length < 8:
            return False, f"Saturated fatty acid chain too short ({chain_length} carbons; require at least 8)"
    
    # 10. Compose a reasoning message.
    reason = (f"Found a terminal carboxylic acid group (acid carbon idx {acid_carbon_idx}) with "
              f"an alpha carbon (idx {alpha_carbon.GetIdx()}) bearing a unique hydroxyl substituent, "
              f"and a predominantly aliphatic chain (chain length = {chain_length} carbons"
              f"{', unsaturated' if unsaturation else ', saturated'}), "
              f"classifying the molecule as a 2–hydroxy fatty acid.")
    return True, reason


# Example usage:
if __name__ == "__main__":
    # Example: (2R)-2-hydroxytetradecanoic acid
    test_smiles = "CCCCCCCCCCCC[C@@H](O)C(O)=O"
    result, classification_reason = is_2_hydroxy_fatty_acid(test_smiles)
    print(result, classification_reason)