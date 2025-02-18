"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: quinone
Definition: Compounds having a fully conjugated cyclic dione structure, such as that of benzoquinones,
derived from aromatic compounds by conversion of an even number of -CH= groups into -C(=O)- groups
(with any necessary rearrangement of double bonds). (Polycyclic and heterocyclic analogues are included.)
"""

from rdkit import Chem

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    
    The strategy is as follows:
      1. Parse the SMILES.
      2. Iterate over each ring in the molecule.
      3. For each ring:
            • Require that the ring is “conjugated” (most atoms have sp2 character).
            • Count candidate carbonyl centers: ring member carbons that are double-bonded to an oxygen.
            • Compute the fraction of carbon atoms in the ring. (Genuine quinones (like benzoquinone) 
              are derived from an all‐carbon aromatic ring.)
      4. If a ring has at least 2 carbonyl groups and most (>=80%) of its atoms are carbon, 
         we classify the molecule as a quinone.
    
    This extra check (the carbon fraction) helps weed out many heterocyclic diones (e.g. uracil‐like rings)
    that are not true quinones.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a quinone, otherwise False.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Loop over each ring
    for ring in rings:
        # We consider rings with at least 5 atoms (most conjugated rings are 5- or 6-membered; quinones are typically 6-membered)
        if len(ring) < 5:
            continue
        
        # Check that most atoms seem sp2-ish.
        sp2_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Many quinone-type rings are either aromatic or have trigonal planar (sp2) geometry.
            if atom.GetIsAromatic() or atom.GetHybridization().name == "SP2":
                sp2_count += 1
        # Require at least 80% of the atoms in the ring be aromatic or sp2.
        if sp2_count / len(ring) < 0.8:
            continue

        # Count how many ring atoms are carbonyl centers.
        carbonyl_count = 0
        carbon_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                carbon_count += 1
                # Look at bonds from this atom
                for bond in atom.GetBonds():
                    # Consider only double bonds.
                    if bond.GetBondType() != Chem.BondType.DOUBLE:
                        continue
                    neighbor = bond.GetOtherAtom(atom)
                    if neighbor.GetAtomicNum() == 8:
                        # Candidate carbonyl: the atom is in the ring, the oxygen typically is exocyclic.
                        carbonyl_count += 1
                        break  # Count once per carbonyl center

        # Check ratio of carbon atoms in ring. (E.g. benzoquinones are derived from all‐carbon rings.)
        if len(ring) > 0:
            carbon_ratio = carbon_count / len(ring)
        else:
            carbon_ratio = 0

        if carbonyl_count >= 2 and carbon_ratio >= 0.8:
            return True, (f"Found ring with atoms {ring}: {carbonyl_count} carbonyl groups and "
                          f"{carbon_count}/{len(ring)} carbon atoms (ratio {carbon_ratio:.2f}).")
    
    return False, "No fully conjugated ring with at least two carbonyl groups and high carbon content was found."


# Example calls (uncomment to test)
# print(is_quinone("Oc1ccc2C(=O)c3c(O)ccc(O)c3C(=O)c2c1O"))  # quinalizarin: a true quinone
# print(is_quinone("COC1=C(C=C2C(=O)C(C)=C(C\\C=C(/C)CC\\C=C(/C)CC\\C=C(/C)CC)C2=O)C1=O"))  # An example structurally similar to quinones