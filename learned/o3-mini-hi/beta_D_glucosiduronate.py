"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate

Definition:
  A carbohydrate acid derivative anion obtained by deprotonation of the carboxyl group
  of any beta-D-glucuronic acid; major species at pH 7.3.

Our revised algorithm requires that the molecule contain a single six-membered (pyranose-like) ring 
composed of one oxygen and five sp3 carbons, with exactly five chiral ring carbons and exactly three 
hydroxyl (-OH) substituents on those carbons. In addition, at least one ring carbon must carry an exocyclic 
carboxylate group (C(=O)[O-]), corresponding to oxidation of the primary alcohol (C6) in D-glucuronic acid.
"""

from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if the molecule contains a beta-D-glucosiduronate fragment.
    
    Checks for:
      - A valid molecule.
      - At least one six-membered ring having exactly one oxygen and five sp3 carbons.
      - Exactly five chiral centers among the ring carbons.
      - Exactly three hydroxyl (-OH) substituents attached to ring carbons.
      - At least one exocyclic carboxylate group (C(=O)[O-]) attached to one of these ring carbons.
    
    Args:
        smiles (str): The SMILES string of the molecule.
    
    Returns:
        bool: True if the fragment is found, False otherwise.
        str: Explanation for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings in the molecule"
    
    # Pre-calculate chiral centers for the entire molecule.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_atoms = {idx for idx, tag in chiral_centers}

    # Helper function to determine if an atom (given its context) is a carboxylate carbon.
    def is_carboxylate_carbon(carbon_atom):
        if carbon_atom.GetAtomicNum() != 6:
            return False
        has_dblO = False
        has_singleO_neg = False
        for nbr in carbon_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
                # Check for a double bond (C=O)
                if bond is not None and bond.GetBondTypeAsDouble() == 2:
                    has_dblO = True
                # And a single-bonded oxygen with -1 charge (O-)
                elif bond is not None and bond.GetBondType() == Chem.BondType.SINGLE:
                    if nbr.GetFormalCharge() == -1:
                        has_singleO_neg = True
        return has_dblO and has_singleO_neg

    # Loop over each ring
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # Only consider six-membered rings.
        
        # Identify ring oxygen and carbons.
        ring_oxygens = []
        ring_carbons = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxygens.append(idx)
            elif atom.GetAtomicNum() == 6:
                ring_carbons.append(idx)
        
        # For a pyranose ring, expect exactly one oxygen and five carbons.
        if len(ring_oxygens) != 1 or len(ring_carbons) != 5:
            continue
        
        # Ensure the ring carbons are sp3 and non‐aromatic.
        valid_carbons = True
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() != Chem.HybridizationType.SP3 or atom.GetIsAromatic():
                valid_carbons = False
                break
        if not valid_carbons:
            continue
        
        # Check chiral centers: require exactly five chiral centers on the ring carbons.
        ring_chiral_count = sum(1 for idx in ring_carbons if idx in chiral_atoms)
        if ring_chiral_count != 5:
            continue
        
        # Count free hydroxyl substituents on ring carbons.
        oh_count = 0
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Skip atoms that are part of the ring.
                if nbr.GetIdx() in ring:
                    continue
                # Look for oxygen substituents that carry at least one hydrogen.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
                    oh_count += 1
        if oh_count != 3:
            continue
        
        # Check that at least one ring carbon is attached to an exocyclic carboxylate.
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # Only consider substituents outside the ring.
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 6 and is_carboxylate_carbon(nbr):
                    return True, ("Found a sugar-like six-membered (pyranose) ring with exactly 5 chiral centers, "
                                  "3 hydroxyl groups, and an exocyclic carboxylate group consistent with beta-D-glucosiduronate")
    
    return False, "No beta-D-glucosiduronate fragment found"

# Example usage:
if __name__ == "__main__":
    # Minimal beta-D-glucuronate (deprotonated) fragment.
    example_smiles = "O[C@@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]"
    result, reason = is_beta_D_glucosiduronate(example_smiles)
    print("Result:", result)
    print("Reason:", reason)