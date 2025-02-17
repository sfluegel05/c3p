"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: triradylglycerol (a glycerol derivative having one substituent group – either acyl, alkyl,
or alk-1-enyl – on each of the three positions sn-1, sn-2, or sn-3).
Functional parent: glycerol (CHEBI:17754), and a type of glycerolipid (CHEBI:35741).
Note that triglycerides (CHEBI:17855) are a subset of triradylglycerols.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES.
    It checks for a simple glycerol backbone (CH2-CH-CH2) where each carbon has exactly one oxygen substituent.
    Then it verifies that each oxygen substituent is connected to a heavy atom (carbon) that appears as one
    of the allowed groups:
      - acyl (the carbon shows at least one double-bond to an oxygen, a common ester pattern),
      - alkyl (a saturated sp3 carbon), or
      - alk-1-enyl (an sp2 carbon, part of a C=C double bond).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as a triradylglycerol, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a simple glycerol backbone (CH2-CH-CH2).
    # This is a simplified approach and may not capture all stereochemical representations.
    glycerol_smarts = "[CH2]-[CH]-[CH2]"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    backbone_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not backbone_matches:
        return False, "No simple glycerol backbone (CH2-CH-CH2) found"
    
    # Initialize an overall list to capture reasons from each candidate backbone match.
    overall_reasons = []
    
    # Iterate over each candidate glycerol backbone
    for match in backbone_matches:
        valid_backbone = True
        candidate_reasons = []  # local reasons for this backbone candidate
        oxygens_found = []
        # Each backbone carbon should have exactly one oxygen substituent (not counting backbone atoms)
        for atom_idx in match:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Count oxygen neighbors that are not part of the backbone
            oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in match and nbr.GetAtomicNum() == 8]
            if len(oxy_neighbors) != 1:
                valid_backbone = False
                candidate_reasons.append(f"Atom idx {atom_idx} has {len(oxy_neighbors)} oxygen substituents (expected 1)")
                break
            oxygens_found.append(oxy_neighbors[0])
            
        if not valid_backbone:
            overall_reasons.append("Candidate backbone failed substituent oxygen check: " + "; ".join(candidate_reasons))
            continue  # try next backbone candidate
        
        # Now check the group attached via each oxygen.
        valid_substituents = True
        for oxy in oxygens_found:
            # Find the substituent heavy atom attached to the oxygen that is not the backbone carbon.
            attached_atoms = [nbr for nbr in oxy.GetNeighbors() if nbr.GetAtomicNum() != 1]  # ignore hydrogens
            substituent = None
            for nbr in attached_atoms:
                # if neighbor is not an oxygen of a glycerol backbone (should be the group fragment)
                if nbr.GetAtomicNum() == 6:
                    substituent = nbr
                    break
            if substituent is None:
                valid_substituents = False
                candidate_reasons.append("Oxygen substituent without a proper carbon neighbor")
                break
            
            # Determine the type of substituent.
            is_acyl = False
            # Check for acyl group: look if there is at least one double bond (C=O) from the substituent carbon.
            for bond in substituent.GetBonds():
                if bond.GetBondType() == rdchem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(substituent)
                    if other.GetAtomicNum() == 8:
                        is_acyl = True
                        break
            
            # Get hybridization to guess alkyl or alk-1-enyl.
            hyb = substituent.GetHybridization()
            
            # Acceptable conditions: either acyl, alkyl (sp3) or alk-1-enyl (sp2)
            if not (is_acyl or (hyb == rdchem.HybridizationType.SP3) or (hyb == rdchem.HybridizationType.SP2)):
                valid_substituents = False
                candidate_reasons.append("Substituent attached via oxygen is not recognized as acyl, alkyl, or alk-1-enyl")
                break
                
        if valid_substituents:
            return True, "Glycerol backbone with three substituents (acyl, alkyl, or alk-1-enyl) recognized"
        else:
            overall_reasons.append("Candidate backbone substituent check failed: " + "; ".join(candidate_reasons))
    
    if overall_reasons:
        return False, " ; ".join(overall_reasons)
    else:
        return False, "No valid glycerol backbone with three suitable substituents found"


# (Optional) Example usage:
if __name__ == '__main__':
    # Test one of the provided TG examples
    test_smiles = "O([C@@H](COC(=O)CCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)C(=O)CCCCCCCCCCCCCC"
    result, reason = is_triradylglycerol(test_smiles)
    print(result, reason)