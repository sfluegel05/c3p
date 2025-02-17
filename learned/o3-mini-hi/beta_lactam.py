"""
Classifies: CHEBI:35627 beta-lactam
"""
"""
Classifies: beta-lactam
A beta-lactam is defined as a lactam in which the amide bond is contained
within a four-membered ring that includes the amide nitrogen and the carbonyl carbon.
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    
    The beta-lactam structure is defined as a 4-membered ring containing an amide bond:
    it must contain a nitrogen atom and a carbonyl carbon (a carbon double-bonded to oxygen)
    with a direct bond between them within the ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is a beta-lactam, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo().AtomRings()
    
    # Loop through rings to check for a four-membered ring with the beta-lactam motif
    for ring in ring_info:
        if len(ring) != 4:
            continue  # We only consider 4-membered rings
        
        # Identify if the ring contains at least one nitrogen
        nitrogen_indices = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7]
        if not nitrogen_indices:
            continue  # No nitrogen in this ring
        
        # Identify carbon atoms in the ring that might be a carbonyl carbon.
        # A carbonyl carbon must be an sp2 carbon that is double-bonded to an oxygen.
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Check that the atom is carbon
            if atom.GetAtomicNum() != 6:
                continue
            # Look for a double bond to oxygen
            found_carbonyl = False
            for bond in atom.GetBonds():
                # Check if bond is a double bond and to an oxygen
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    other = bond.GetOtherAtom(atom)
                    if other.GetAtomicNum() == 8:
                        found_carbonyl = True
                        break
            if not found_carbonyl:
                continue  # This carbon is not carbonyl
            
            # At this point, 'atom' is a candidate for the carbonyl carbon.
            # Check if it is directly bonded to a nitrogen within the same ring.
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetAtomicNum() == 7 and neighbor.GetIdx() in ring:
                    # Found a nitrogen directly bonded to the carbonyl carbon in the ring.
                    return True, "Beta-lactam ring detected: 4-membered ring with a carbonyl carbon and an adjacent nitrogen"
    
    return False, "No beta-lactam ring found (4-membered ring with the required amide bond was not detected)"

# Example usage:
if __name__ == "__main__":
    # Example beta-lactam: azetidin-2-one is the simplest beta-lactam
    test_smiles = "O=C1CCN1"
    result, reason = is_beta_lactam(test_smiles)
    print("Result:", result)
    print("Reason:", reason)