"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate

Definition:
  A carbohydrate acid derivative anion obtained by deprotonation of the carboxyl group
  of any beta-D-glucuronic acid; major species at pH 7.3.

Our revised algorithm searches for a pyranose-like ring that meets all the criteria:
  - A six-membered ring containing exactly one sp³ oxygen and five sp³ (nonaromatic) carbons.
  - Exactly five chiral centers among the ring carbons.
  - Exactly three free hydroxyl (-OH) substituents attached to those ring carbons.
  - Exactly one exocyclic carboxylate group (pattern C(=O)[O-]) attached to one of the ring carbons.
"""

from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Check if the molecule contains a beta-D-glucosiduronate fragment.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a beta-D-glucosiduronate fragment is found, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings in the molecule"

    # Pre-calculate chiral centers over the entire molecule.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_atoms = {idx for idx, _ in chiral_centers}

    # Helper to check if an atom is a carboxylate carbon (C(=O)[O-])
    def is_carboxylate_carbon(carbon_atom):
        if carbon_atom.GetAtomicNum() != 6:
            return False
        has_dblO = False
        has_singleO_neg = False
        for nbr in carbon_atom.GetNeighbors():
            # We expect one double-bonded oxygen and one single-bonded oxygen with -1 charge.
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                # For a double-bonded oxygen
                if bond.GetBondTypeAsDouble() == 2:
                    has_dblO = True
                # For a single-bonded oxygen that carries a -1 formal charge.
                elif bond.GetBondType() == Chem.BondType.SINGLE and nbr.GetFormalCharge() == -1:
                    has_singleO_neg = True
        return has_dblO and has_singleO_neg

    # Loop over each ring in the molecule.
    for ring in ring_info.AtomRings():
        # Only consider six-membered rings.
        if len(ring) != 6:
            continue

        # Partition ring atoms into oxygens and carbons.
        ring_oxygens = []
        ring_carbons = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxygens.append(idx)
            elif atom.GetAtomicNum() == 6:
                ring_carbons.append(idx)

        # For a pyranose ring candidate, we require exactly one oxygen and five carbons.
        if len(ring_oxygens) != 1 or len(ring_carbons) != 5:
            continue

        # Check that the ring oxygen is sp3.
        ring_oxy_atom = mol.GetAtomWithIdx(ring_oxygens[0])
        if ring_oxy_atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue

        # Check that all ring carbons are sp3 and nonaromatic.
        valid_carbons = True
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() != Chem.HybridizationType.SP3 or atom.GetIsAromatic():
                valid_carbons = False
                break
        if not valid_carbons:
            continue

        # Count chiral centers on the ring carbons; require exactly five.
        ring_chiral_count = sum(1 for idx in ring_carbons if idx in chiral_atoms)
        if ring_chiral_count != 5:
            continue

        # Count free hydroxyl (-OH) substituents on ring carbons.
        oh_count = 0
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Ignore connections within the ring.
                if nbr.GetIdx() in ring:
                    continue
                # Look for oxygen substituents that have at least one hydrogen.
                # (Carboxylate oxygens usually have no hydrogen.)
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
                    oh_count += 1
        if oh_count != 3:
            continue

        # Count the exocyclic carboxylate groups attached to ring carbons.
        carboxylate_count = 0
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Process only substituents not in the ring.
                if nbr.GetIdx() in ring:
                    continue
                # Look for a carbon that might be a carboxylate.
                if nbr.GetAtomicNum() == 6 and is_carboxylate_carbon(nbr):
                    carboxylate_count += 1
        if carboxylate_count != 1:
            continue

        # If all checks pass, we have found a beta-D-glucosiduronate fragment.
        return True, ("Found a pyranose-like ring with 5 sp3 chiral carbons, 3 hydroxyl substituents, "
                      "and exactly 1 exocyclic carboxylate group consistent with beta-D-glucosiduronate")

    return False, "No beta-D-glucosiduronate fragment found"

# Example usage:
if __name__ == "__main__":
    # Minimal deprotonated beta-D-glucuronate fragment example.
    example_smiles = "O[C@@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]"
    result, reason = is_beta_D_glucosiduronate(example_smiles)
    print("Result:", result)
    print("Reason:", reason)