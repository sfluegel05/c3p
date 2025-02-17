"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate

Definition:
  A carbohydrate acid derivative anion obtained by deprotonation of the carboxy group
  of any beta-D-glucuronic acid; major species at pH 7.3.

This implementation first looks for candidate six-membered (pyranose-like) rings
composed of five sp3 carbons and one oxygen. It then requires at least three chiral
centers and three –OH substituents on the ring carbons (features typical of a sugar)
plus the attachment of a carboxylate group (C(=O)[O-]) off one of the ring carbons.
"""

from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if the molecule has a beta-D-glucosiduronate fragment.
    
    The simplified algorithm checks for:
      - A valid molecule.
      - At least one six-membered ring that appears sugar-like, meaning:
          * It has exactly six atoms with one oxygen and five carbons.
          * All ring carbons are sp3 and (usually) chiral – we require at least 3 chiral centers.
          * At least three hydroxyl (-OH) substituents attached to ring carbons.
      - An exocyclic carboxylate group (C(=O)[O-]) attached to one of the ring carbons.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a beta-D-glucosiduronate-like fragment is found, False otherwise.
        str: Reason for the classification.
    """
    # Parse SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring info.
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings in the molecule"
    
    # Pre-calculate chiral centers for the whole molecule.
    # (Returns list of (atom index, chirality tag))
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_atoms = {idx for idx, tag in chiral_centers}
    
    # Define a helper to decide if a given carbon atom (with mol context) is part of a carboxylate group.
    def is_carboxylate_carbon(carbon_atom):
        # Must be carbon.
        if carbon_atom.GetAtomicNum() != 6:
            return False
        dbl_o = False
        single_o = False
        for nbr in carbon_atom.GetNeighbors():
            # We require one double-bonded oxygen and one single bond to oxygen bearing a -1 charge.
            if nbr.GetAtomicNum() == 8:
                bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
                if bond is None:
                    continue
                if bond.GetBondTypeAsDouble() == 2:
                    dbl_o = True
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    if nbr.GetFormalCharge() == -1:
                        single_o = True
        return dbl_o and single_o

    # Loop over each ring in the molecule
    for ring in ring_info.AtomRings():
        if len(ring) != 6:
            continue  # only consider six-membered rings
        
        # Count ring oxygens and collect ring carbon indices.
        ring_oxygens = []
        ring_carbons = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxygens.append(idx)
            elif atom.GetAtomicNum() == 6:
                ring_carbons.append(idx)
        
        # For a pyranose ring we expect exactly one oxygen and five carbons.
        if len(ring_oxygens) != 1 or len(ring_carbons) != 5:
            continue
        
        # Check that the ring carbons are sp3, non‐aromatic.
        sp3_carbons = True
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() != Chem.HybridizationType.SP3 or atom.GetIsAromatic():
                sp3_carbons = False
                break
        if not sp3_carbons:
            continue
        
        # Count how many of these ring carbons are chiral.
        chiral_count = sum(1 for idx in ring_carbons if idx in chiral_atoms)
        if chiral_count < 3:
            continue  # not sugar-like enough
        
        # Count hydroxyl substituents on ring carbons.
        oh_count = 0
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            # Look for oxygen neighbors (not in ring) that are -OH.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # check that this oxygen is attached as -OH (has at least one hydrogen)
                    # We check total number of hydrogens (implicit+explicit)
                    if nbr.GetTotalNumHs() >= 1:
                        oh_count += 1
        if oh_count < 3:
            continue  # too few hydroxyl substituents for a typical sugar
        
        # Finally, check that at least one ring carbon has an exocyclic carboxylate group.
        for idx in ring_carbons:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # The candidate must be a carbon that is a carboxylate.
                if nbr.GetAtomicNum() == 6 and is_carboxylate_carbon(nbr):
                    return True, ("Found a sugar-like six-membered (pyranose) ring with %d chiral centers "
                                  "and %d hydroxyl groups and an exocyclic carboxylate group "
                                  "consistent with beta-D-glucosiduronate" % (chiral_count, oh_count))
    
    return False, "No beta-D-glucosiduronate fragment found"

# Example usage:
if __name__ == "__main__":
    # This SMILES is for a minimal beta-D-glucuronate fragment.
    example_smiles = "O[C@@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]"
    result, reason = is_beta_D_glucosiduronate(example_smiles)
    print("Result:", result)
    print("Reason:", reason)