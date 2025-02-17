"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: triradylglycerol – a glycerol derivative having one substituent group 
(at sn-1, sn-2, and sn-3) that is either acyl, alkyl, or alk-1-enyl.
Functional parent: glycerol (CHEBI:17754);
Child: triglyceride (CHEBI:17855); Parent: glycerolipid (CHEBI:35741).

Note: This classifier first eliminates molecules with phosphorus (such as common glycerophospholipids)
and then looks for a simple glycerol backbone (CH2-CH-CH2) where every carbon has exactly one oxygen 
substituent. The oxygen must be connected to a carbon (and not to phosphorus or other atoms). Finally, 
the attached carbon is examined to decide if it is acyl (presence of a C=O bond) or, if not, whether 
its hybridization is sp3 (alkyl) or sp2 (alk-1-enyl).
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES.
    
    The algorithm first rejects molecules that contain phosphorus (P), then searches for a glycerol 
    backbone represented by a simple SMARTS pattern [CH2]-[CH]-[CH2]. For any matching backbone, each 
    of the three carbons is required to have exactly one oxygen substituent that is not part of the backbone 
    and that oxygen must then be connected to a carbon (and not to any disallowed element). The carbon attached 
    via oxygen is then evaluated:
      - If it has at least one double bond to oxygen (C=O) it is considered an acyl group.
      - Otherwise its hybridization must be either sp3 (alkyl) or sp2 (alk-1-enyl) in order to be accepted.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a triradylglycerol, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules that contain phosphorus (or other non-carbon heteroatoms in substituents if needed).
    if any(atom.GetAtomicNum() == 15 for atom in mol.GetAtoms()):
        return False, "Contains phosphorus (likely a glycerophospholipid), not a triradylglycerol"
    
    # Define a SMARTS pattern for a simple glycerol backbone (CH2-CH-CH2)
    glycerol_smarts = "[CH2]-[CH]-[CH2]"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    backbone_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not backbone_matches:
        return False, "No simple glycerol backbone (CH2-CH-CH2) found"
    
    # Try each candidate glycerol backbone match
    overall_reasons = []
    for match in backbone_matches:
        valid_backbone = True
        candidate_reasons = []
        substituent_atoms = []  # will hold the carbon attached to the oxygen in each case
        
        # For each atom in the candidate backbone, check that it has exactly one oxygen neighbor (not part of the backbone)
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in match and nbr.GetAtomicNum() == 8]
            if len(oxy_neighbors) != 1:
                valid_backbone = False
                candidate_reasons.append(f"Backbone atom idx {atom_idx} has {len(oxy_neighbors)} oxygen substituents (expected 1)")
                break
            oxy_atom = oxy_neighbors[0]
            # Check that oxygen is not connected to phosphorus or other disallowed elements
            attached_heavy = [nbr for nbr in oxy_atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != atom.GetIdx()]
            if not attached_heavy:
                valid_backbone = False
                candidate_reasons.append(f"Oxygen at atom idx {oxy_atom.GetIdx()} is not attached to any heavy atom besides backbone carbon")
                break
            # We expect the substituent to start with a carbon
            substituent = None
            for nbr in attached_heavy:
                if nbr.GetAtomicNum() == 6:
                    substituent = nbr
                    break
            if substituent is None:
                valid_backbone = False
                candidate_reasons.append("Substituent group attached via oxygen does not start with carbon")
                break
            substituent_atoms.append(substituent)
        
        if not valid_backbone:
            overall_reasons.append("Candidate backbone failed: " + "; ".join(candidate_reasons))
            continue  # try next candidate backbone
        
        # Now check for each substituent group if it qualifies as acyl, alkyl, or alk-1-enyl.
        for sub_atom in substituent_atoms:
            qualifies = False
            # Check for acyl: look for at least one double bond (C=O) beginning at this carbon.
            for bond in sub_atom.GetBonds():
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(sub_atom)
                    if other.GetAtomicNum() == 8:
                        qualifies = True  # found a C=O bond → acyl
                        break
            if not qualifies:
                # If not acyl, then check the hybridization:
                hyb = sub_atom.GetHybridization()
                if hyb == rdchem.HybridizationType.SP3 or hyb == rdchem.HybridizationType.SP2:
                    qualifies = True
            if not qualifies:
                valid_backbone = False
                candidate_reasons.append("A substituent does not appear to be acyl, alkyl, or alk-1-enyl")
                break
        
        if valid_backbone:
            return True, "Glycerol backbone with three substituents (acyl, alkyl, or alk-1-enyl) recognized"
        else:
            overall_reasons.append("Candidate backbone substituent check failed: " + "; ".join(candidate_reasons))
    
    return False, " ; ".join(overall_reasons)

# (Optional) Example usage:
if __name__ == '__main__':
    # Test one TG example (this one should be classified as a triradylglycerol)
    test_smiles = "O([C@@H](COC(=O)CCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)C(=O)CCCCCCCCCCCCCC"
    result, reason = is_triradylglycerol(test_smiles)
    print(result, reason)