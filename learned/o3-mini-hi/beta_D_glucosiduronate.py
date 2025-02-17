"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: beta-D-glucosiduronate

Definition:
  A carbohydrate acid derivative anion obtained by deprotonation of the carboxy group 
  of any beta-D-glucuronic acid; major species at pH 7.3.

Revised algorithm:
  1. Look for six-membered rings that contain exactly one sp3 oxygen and five sp3 nonaromatic carbons.
  2. Check that each of the five ring carbons is chiral.
  3. For each ring carbon, require that it has exactly one exocyclic neighbor (i.e. not part of the ring) 
     that is either an oxygen atom (expected for free hydroxyl or glycosidic link) or a carbon (if it is a carboxylate).
  4. For each such substituent:
     • Classify as "OH" if it is oxygen and has at least one hydrogen.
     • Classify as "Gly" if it is oxygen and carries no hydrogens.
     • Classify as "COO" if it is a carbon that fits the carboxylate pattern: bonded via a double bond 
       to an oxygen and via a single bond to an oxygen bearing a –1 formal charge.
  5. Accept the candidate ring only if exactly one ring carbon shows a carboxylate substituent,
     exactly three show a free –OH and exactly one shows a glycosidic substituent (the anomeric center).
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
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No rings in the molecule"
    
    # Pre-calculate chiral centers on entire molecule.
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    chiral_atoms = {idx for idx, _ in chiral_centers}

    # Helper function: Test if an atom is a carboxylate carbon.
    def is_carboxylate_carbon(carbon_atom):
        # It must be a carbon.
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
            # Check for double bond to oxygen.
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                dblO_found = True
            # Check single bond to oxygen with formal charge -1.
            elif bond.GetBondType() == Chem.BondType.SINGLE and nbr.GetFormalCharge() == -1:
                singleO_neg_found = True
        return dblO_found and singleO_neg_found

    # Loop over candidate rings
    for ring in ring_info.AtomRings():
        # Only consider rings of size 6.
        if len(ring) != 6:
            continue

        # Partition atoms in the ring into oxygen(s) and carbons.
        ring_oxy_indices = []
        ring_carbon_indices = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 8:
                ring_oxy_indices.append(idx)
            elif atom.GetAtomicNum() == 6:
                ring_carbon_indices.append(idx)
        # For a pyranose candidate, expect exactly 1 ring oxygen and 5 ring carbons.
        if len(ring_oxy_indices) != 1 or len(ring_carbon_indices) != 5:
            continue
        
        # Check that the ring oxygen is sp3.
        ring_oxy = mol.GetAtomWithIdx(ring_oxy_indices[0])
        if ring_oxy.GetHybridization() != Chem.HybridizationType.SP3:
            continue

        # Check that all ring carbons are sp3 and nonaromatic.
        valid_carbons = True
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetHybridization() != Chem.HybridizationType.SP3 or atom.GetIsAromatic():
                valid_carbons = False
                break
        if not valid_carbons:
            continue

        # Require that each of the five ring carbons is chiral.
        ring_chiral = sum(1 for idx in ring_carbon_indices if idx in chiral_atoms)
        if ring_chiral != 5:
            continue

        # For each ring carbon, collect its exocyclic substituents (neighbors that are not in the ring)
        substituent_types = {}  # key: ring carbon index, value: type "OH", "Gly" or "COO"
        valid_candidate = True
        for idx in ring_carbon_indices:
            atom = mol.GetAtomWithIdx(idx)
            ext_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in ring]
            # Each ring carbon should have exactly one exocyclic substituent.
            if len(ext_neighbors) != 1:
                valid_candidate = False
                break
            
            nbr = ext_neighbors[0]
            # If neighbor is oxygen
            if nbr.GetAtomicNum() == 8:
                # Count explicit hydrogens (the sum of implicit and explicit)
                num_h = nbr.GetTotalNumHs()
                # If this oxygen carries one or more hydrogen, we consider it a free hydroxyl.
                if num_h >= 1:
                    substituent_types[idx] = "OH"
                else:
                    # No hydrogen means likely a glycosidic (substituted) oxygen.
                    substituent_types[idx] = "Gly"
            # Otherwise if neighbor is carbon, check if it qualifies as a carboxylate.
            elif nbr.GetAtomicNum() == 6:
                if is_carboxylate_carbon(nbr):
                    substituent_types[idx] = "COO"
                else:
                    # Not a carboxylate; not expected.
                    valid_candidate = False
                    break
            else:
                valid_candidate = False
                break
        if not valid_candidate:
            continue

        # Now check that among the five ring carbons the substituents count as required:
        # exactly one "COO", three "OH" and one "Gly".
        oh_count  = sum(1 for t in substituent_types.values() if t == "OH")
        coo_count = sum(1 for t in substituent_types.values() if t == "COO")
        gly_count = sum(1 for t in substituent_types.values() if t == "Gly")
        if oh_count == 3 and coo_count == 1 and gly_count == 1:
            return True, ("Found a pyranose-like ring with 5 sp³ chiral carbons, "
                          "exactly one carboxylate substituent, three free hydroxyls, "
                          "and one glycosidically substituted (anomeric) center consistent with beta-D-glucosiduronate")
    
    return False, "No beta-D-glucosiduronate fragment found"

# Example usage:
if __name__ == "__main__":
    # Minimal beta-D-glucuronide example (already deprotonated)
    example_smiles = "O[C@@H]1O[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]"
    result, reason = is_beta_D_glucosiduronate(example_smiles)
    print("Result:", result)
    print("Reason:", reason)