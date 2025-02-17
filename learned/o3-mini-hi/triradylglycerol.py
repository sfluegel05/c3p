"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: triradylglycerol, defined as a glycerol compound having one of three possible substituent groups – 
either acyl, alkyl, or alk-1-enyl – at each of the three positions (sn-1, sn-2, sn-3). Its functional parent is glycerol 
(CHEBI:17754) and it is a type of glycerolipid (CHEBI:35741). Note that triglycerides (CHEBI:17855) are a subset of this class.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol has a glycerol backbone (i.e. a CH2-CH-CH2 chain) with each of the three carbons
    substituted by an oxygen that is linked to either an acyl (ester), alkyl (ether) or alk-1-enyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a triradylglycerol, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a simple glycerol backbone:
    # three carbons: first and third are CH2, middle is CH.
    glycerol_smarts = "[CH2]-[CH]-[CH2]"
    glycerol_pattern = Chem.MolFromSmarts(glycerol_smarts)
    backbone_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not backbone_matches:
        return False, "No simple glycerol backbone (CH2-CH-CH2) found"
    
    # For each candidate backbone, check if each carbon has exactly one oxygen substituent outside the chain.
    for match in backbone_matches:
        # match gives indices of the three carbons in order.
        oxygens_found = []
        valid_backbone = True
        for i, atom_idx in enumerate(match):
            atom = mol.GetAtomWithIdx(atom_idx)
            oxygens = []
            # Examine neighbors of the glycerol carbon
            for nbr in atom.GetNeighbors():
                # Skip if neighbor is one of the backbone carbons
                if nbr.GetIdx() in match:
                    continue
                # We require the substituent to be an oxygen.
                if nbr.GetAtomicNum() == 8:
                    oxygens.append(nbr)
            # For each backbone carbon, we require exactly one oxygen substituent.
            if len(oxygens) != 1:
                valid_backbone = False
                break
            oxygens_found.append(oxygens[0])
        if not valid_backbone:
            continue  # try the next backbone match
        
        # At this point we have a candidate glycerol backbone with three oxygen substituents.
        # Now we check each substituent oxygen to ensure it links to a heavy fragment that is either:
        #   - an acyl group (ester): oxygen attached to a carbon that carries a carbonyl double bond,
        #   - an alkyl group (saturated, sp3 carbon) or
        #   - an alk-1-enyl group (oxygen attached to an sp2 carbon involved in a C=C bond).
        valid_substituents = True
        reason_msgs = []
        for oxy in oxygens_found:
            # Get the neighbor attached to oxygen that is not the glycerol carbon.
            nbrs = [n for n in oxy.GetNeighbors() if n.GetAtomicNum() != 1]  # ignore hydrogens
            # It should be exactly one neighbor aside from the glycerol attachment.
            # (In many cases additional hydrogens are implicit; we use heavy-atom connectivity.)
            substituent = None
            for nbr in nbrs:
                # if this neighbor is not part of our backbone, choose it as substituent
                if nbr.GetIdx() not in match:
                    substituent = nbr
                    break
            if substituent is None:
                valid_substituents = False
                reason_msgs.append("Oxygen substituent without a heavy atom neighbor")
                break

            # Now check type of substituent: acyl (ester), alkyl (saturated sp3) or alk-1-enyl.
            # Acyl group: Look for substituent carbon that has at least one double-bonded oxygen.
            if substituent.GetAtomicNum() == 6:
                is_acyl = False
                is_alkenyl = False
                # Check each neighbor bond of substituent for a double bond to O.
                for bond in mol.GetBonds():
                    # We look for a bond from this substituent to an O with bond type DOUBLE.
                    if (bond.GetBeginAtomIdx() == substituent.GetIdx() or bond.GetEndAtomIdx() == substituent.GetIdx()):
                        other = bond.GetOtherAtom(substituent)
                        if other.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                            is_acyl = True
                            break
                # If not acyl, then check hybridization.
                hyb = substituent.GetHybridization()
                if not is_acyl:
                    # If sp2 then assume alk-1-enyl, if sp3 then assume alkyl.
                    if hyb == rdchem.HybridizationType.SP2:
                        is_alkenyl = True
                # If none of the conditions met, then fail this substituent.
                if not (is_acyl or (hyb == rdchem.HybridizationType.SP3) or is_alkenyl):
                    valid_substituents = False
                    reason_msgs.append("Substituent attached via oxygen does not appear as acyl, alkyl, or alk-1-enyl")
                    break
            else:
                valid_substituents = False
                reason_msgs.append("Substituent attached to oxygen is not a carbon")
                break
        
        if valid_substituents:
            return True, "Glycerol backbone with 3 substituents (acyl, alkyl, or alk-1-enyl) recognized"
    
    # If no backbone match passed the substituent checks
    if reason_msgs:
        return False, "; ".join(reason_msgs)
    else:
        return False, "No valid glycerol backbone with three suitable substituents found"

# (Optional) Example usage:
if __name__ == '__main__':
    # Try one of the provided TG examples
    test_smiles = "O([C@@H](COC(=O)CCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)C(=O)CCCCCCCCCCCCCC"
    result, reason = is_triradylglycerol(test_smiles)
    print(result, reason)