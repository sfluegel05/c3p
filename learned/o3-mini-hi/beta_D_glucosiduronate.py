"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate

Definition:
  A carbohydrate acid derivative anion obtained by deprotonation of the carboxyl group 
  of any beta-D-glucuronic acid; major species at pH 7.3.

Revised algorithm:
  1. Search for a six-membered ring having exactly one sp3 oxygen and exactly five sp3, 
     nonaromatic carbons. Require that five chiral centres are present (one per C).
  2. For each ring carbon, inspect its substituents (neighbors not in the ring):
      • Count as a free –OH if the neighbor is oxygen and carries at least one hydrogen.
      • Count as a carboxylate if it is a carbon atom that has one double‐bonded oxygen and one
        single‐bonded oxygen carrying a –1 formal charge.
      • Count as a “glycosidic” oxygen if the neighbor is oxygen and carries no (free) hydrogen.
  3. Then require that among the five ring carbons exactly:
      • One shows (only) an exocyclic carboxylate group,
      • Three show a free –OH substituent,
      • One (the anomeric carbon) shows a non‐OH oxygen substituent (i.e. glycosidic linkage).
  
If a ring meeting these tight conditions is found, we conclude that the SMILES
contains a beta-D-glucosiduronate fragment.
"""

from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if the molecule contains a beta-D-glucosiduronate fragment.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the beta-D-glucosiduronate fragment is found, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings in the molecule"
    
    # Pre-calculate chiral centers in the entire molecule.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_atoms = {idx for idx, _ in chiral_centers}

    # Helper: Check if a given carbon atom is a carboxylate carbon (C(=O)[O-])
    def is_carboxylate_carbon(carbon_atom):
        if carbon_atom.GetAtomicNum() != 6:
            return False
        dblO_found = False
        singleO_neg_found = False
        for nbr in carbon_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 8:
                continue
            bond = mol.GetBondBetweenAtoms(carbon_atom.GetIdx(), nbr.GetIdx())
            if not bond:
                continue
            # Check for a double bond to oxygen.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                dblO_found = True
            # Check for a single bond to an oxygen carrying a -1 charge.
            elif bond.GetBondType() == Chem.BondType.SINGLE and nbr.GetFormalCharge() == -1:
                singleO_neg_found = True
        return dblO_found and singleO_neg_found

    # Loop over each ring.
    for ring in ring_info.AtomRings():
        # Only look at rings of size 6
        if len(ring) != 6:
            continue

        # Partition ring atoms into oxygens and carbons.
        ring_oxy_indices = []
        ring_carbon_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxy_indices.append(idx)
            elif atom.GetAtomicNum() == 6:
                ring_carbon_indices.append(idx)
        
        # For a pyranose candidate: one ring oxygen and five ring carbons.
        if len(ring_oxy_indices) != 1 or len(ring_carbon_indices) != 5:
            continue
        
        # The ring oxygen must be sp3.
        ring_oxy = mol.GetAtomWithIdx(ring_oxy_indices[0])
        if ring_oxy.GetHybridization() != Chem.HybridizationType.SP3:
            continue
        
        # All ring carbons must be sp3 and nonaromatic.
        valid_carbons = True
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() != Chem.HybridizationType.SP3 or atom.GetIsAromatic():
                valid_carbons = False
                break
        if not valid_carbons:
            continue
        
        # Check that exactly 5 chiral centers occur on the ring carbons.
        ring_chiral = sum(1 for idx in ring_carbon_indices if idx in chiral_atoms)
        if ring_chiral != 5:
            continue

        # For each ring carbon, record its exocyclic substituent type.
        # For atoms in the candidate ring we will count one of:
        #   "OH" : free hydroxyl (exocyclic oxygen carrying >=1 H)
        #   "COO" : carboxylate substituent (exocyclic carbon satisfying is_carboxylate_carbon)
        #   "Gly" : substituted oxygen (exocyclic oxygen that carries no H),
        #           typically at the anomeric position.
        substituent_types = {}  # key: ring carbon index, value: type string
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            sub_type = None
            # Look at all neighbors that are not within the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Check if neighbor is oxygen with at least one hydrogen (free OH)
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() >= 1:
                    # We expect each ring carbon to give at most one free OH.
                    sub_type = "OH"
                    break  # highest priority for free hydroxyl.
                # Otherwise if neighbor is oxygen with no hydrogen, mark it as glycosidic.
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() == 0:
                    sub_type = "Gly"
                    # Do not break; a carbon might have two substituents.
                    # But we expect exactly one substituent per ring carbon.
                # If neighbor is carbon, check if it is a carboxylate.
                if nbr.GetAtomicNum() == 6 and is_carboxylate_carbon(nbr):
                    sub_type = "COO"
                    # Carboxylate takes precedence.
                    break
            if sub_type is not None:
                substituent_types[idx] = sub_type
        # Now, among the five ring carbons, count the substituent types.
        oh_count    = sum(1 for t in substituent_types.values() if t == "OH")
        coo_count   = sum(1 for t in substituent_types.values() if t == "COO")
        gly_count   = sum(1 for t in substituent_types.values() if t == "Gly")

        # Our expectation for a beta-D-glucosiduronate fragment (typically in glucuronide conjugates)
        # is that exactly one ring carbon (C5 in the original sugar) is linked to a carboxylate,
        # three ring carbons carry free hydroxyl groups (C2, C3, C4),
        # and the remaining ring carbon (the anomeric center, C1) is glycosidically substituted.
        if oh_count == 3 and coo_count == 1 and gly_count == 1:
            return True, ("Found a pyranose-like ring with 5 sp3 chiral carbons, 3 free hydroxyl substituents, "
                          "1 carboxylate group, and 1 glycosidically substituted anomeric center consistent with beta-D-glucosiduronate")
    
    return False, "No beta-D-glucosiduronate fragment found"

# Example usage:
if __name__ == "__main__":
    # Minimal deprotonated beta-D-glucuronide fragment example.
    example_smiles = "O[C@@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]"
    result, reason = is_beta_D_glucosiduronate(example_smiles)
    print("Result:", result)
    print("Reason:", reason)