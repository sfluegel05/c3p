"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17beta-hydroxy steroid
Definition: A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
This simplified classifier first checks for a steroid‐like tetracyclic fused ring (a common
representation is three six‐membered rings fused with one five‐membered ring) and then looks
for an aliphatic –OH group attached to a carbon that lies in a five‐membered ring (the D‐ring).
Note: Stereochemical details (alpha vs beta) are difficult to infer solely by SMARTS; however,
our approach assumes that an –OH on the D‐ring of a steroid nucleus is likely the 17β–OH.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    The approach is to:
      1. Parse the SMILES.
      2. Search for a steroid-like nucleus by matching a tetracyclic fused ring pattern.
      3. Look for a hydroxyl (-OH) group attached to an sp3 carbon that is part of a five-membered ring,
         as expected for the D-ring containing C17.
      4. Optionally, check that the carbon bearing the –OH is chiral (to indicate defined stereochemistry).
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        (bool, str): True plus a reason if the molecule appears to be a 17beta-hydroxy steroid,
                     False plus a reason otherwise.
                     If SMILES parsing fails, returns (False, "Invalid SMILES string").
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Step 1: Check for a steroid nucleus ---
    # A common (simplified) way to find a steroid core is to look for a fused ring system
    # that contains three six-membered rings and one five-membered ring.
    # The following SMARTS tries to capture part of this fused system:
    steroid_core_smarts = "C1CC2CCC3C(C2)CCC13"
    steroid_core = Chem.MolFromSmarts(steroid_core_smarts)
    if not mol.HasSubstructMatch(steroid_core):
        return False, "Steroid nucleus not found"

    # --- Step 2: Look for a hydroxyl group (-OH) on a saturated carbon ---
    # We require the pattern: an aliphatic carbon (sp3, not part of a carbonyl) singly bonded to an -OH.
    oh_pattern = Chem.MolFromSmarts("[CX4;!$(C=O)][OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if not oh_matches:
        return False, "No appropriate hydroxyl (-OH) group on an aliphatic carbon found"

    # Get ring information for further analysis.
    ring_info = mol.GetRingInfo()
    
    # --- Step 3: Look for a candidate C17 hydroxyl ---
    # In many steroids, the D-ring is five-membered. Thus, we check if the carbon bearing the -OH
    # (the first atom in the pattern match) belongs to a 5-membered ring.
    candidate_found = False
    reason = ""
    for match in oh_matches:
        carbon_idx = match[0]  # the sp3 carbon attached to the OH
        # Check if this carbon is in any 5-membered ring.
        in_five_member_ring = False
        for ring in ring_info.AtomRings():
            if len(ring) == 5 and carbon_idx in ring:
                in_five_member_ring = True
                break
        if in_five_member_ring:
            # Optionally, check that the carbon has a defined chirality.
            atom = mol.GetAtomWithIdx(carbon_idx)
            # If the chiral tag is set (i.e. not unspecified) we assume configuration is defined.
            if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                candidate_found = True
                reason = ("Steroid nucleus found and a hydroxyl group is attached at a chiral carbon "
                          "in a 5-membered ring (likely the D-ring, position C17) with defined stereochemistry.")
                break
            else:
                # Even if chirality is not set, we mark as candidate (with a note).
                candidate_found = True
                reason = ("Steroid nucleus found and a hydroxyl group is attached on a carbon in a 5-membered ring "
                          "(possibly C17), but stereochemistry is not explicitly defined.")
                break

    if not candidate_found:
        return False, "No hydroxyl group found on a carbon in a 5-membered ring (expected for 17β–OH)"
    
    return True, reason