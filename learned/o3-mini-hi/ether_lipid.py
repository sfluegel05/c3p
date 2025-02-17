"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: Ether Lipid
Definition: A lipid similar in structure to a glycerolipid but in which one or more of the carbon atoms on glycerol 
is bonded to an alkyl chain via an ether linkage (C-O-C) rather than the usual ester linkage (O-C(=O)-).
"""

from rdkit import Chem

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    The molecule must contain a glycerol backbone and at least one ether linkage (C-O-C)
    that is not part of an ester group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as an ether lipid, False otherwise.
        str: Explanation for the classification decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # We first search for a glycerol backbone.
    # This simple SMARTS pattern is intended to capture the HO-CH2-CHOH-CH2-OH motif.
    # Note that many lipids contain more decorations (phosphate groups etc.), but this serves as a hint.
    glycerol_smarts = "OCC(O)CO"  # simple representation of glycerol motif
    glycerol_query = Chem.MolFromSmarts(glycerol_smarts)
    glycerol_matches = mol.GetSubstructMatches(glycerol_query)
    if not glycerol_matches:
        return False, "Glycerol backbone not found"

    # For later use, collect all atom indices that are part of any glycerol backbone hit.
    glycerol_atoms = set()
    for match in glycerol_matches:
        for idx in match:
            glycerol_atoms.add(idx)

    # Define a SMARTS pattern for an ether linkage: a carbon-oxygen-carbon unit.
    # This will match any C-O-C bond.
    ether_smarts = "[#6]-O-[#6]"
    ether_query = Chem.MolFromSmarts(ether_smarts)
    ether_matches = mol.GetSubstructMatches(ether_query)
    if not ether_matches:
        return False, "No C-O-C ether linkage found"

    # Function to check if a given carbon atom is part of a carbonyl group.
    def is_carbonyl(carbon):
        # A carbon is considered part of a carbonyl if it has a double bond to an oxygen.
        for bond in carbon.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(carbon)
                if other.GetAtomicNum() == 8:  # oxygen
                    return True
        return False

    # Now loop over found ether linkages and check:
    # (1) They are truly ethers (i.e. not part of an ester)
    # (2) At least one of the carbons in the linkage is part of the glycerol backbone.
    for match in ether_matches:
        # match returns a tuple (c1, o, c2)
        c1 = mol.GetAtomWithIdx(match[0])
        o_atom = mol.GetAtomWithIdx(match[1])
        c2 = mol.GetAtomWithIdx(match[2])
        # Exclude ether linkages that are actually ester bonds.
        # In an ester, one of the two carbons in the C-O-C motif will be a carbonyl carbon.
        if is_carbonyl(c1) or is_carbonyl(c2):
            continue

        # Check if either carbon is in the glycerol backbone.
        if match[0] in glycerol_atoms or match[2] in glycerol_atoms:
            return True, "Molecule contains a glycerol backbone and at least one ether linkage not part of an ester"
    
    # If no ether linkage directly attached to the glycerol backbone is found:
    return False, "No ether linkage found that is attached to the glycerol backbone"

# Example usage:
if __name__ == "__main__":
    # Example SMILES for an ether lipid (one of the examples provided):
    smiles_examples = [
        "P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)COCCCCCCCCCCCCCCCCCC)(O)(O)=O", # PA(O-18:0/22:1(11Z)) has an ether linkage (look for COC... from glycerol)
        "CCCCCCCC\\C=C/CCCCCCCCOC[C@@H](O)COP([O-])(=O)OCC[N+](C)(C)C",  # 1-oleyl-sn-glycero-3-phosphocholine
    ]
    for s in smiles_examples:
        result, reason = is_ether_lipid(s)
        print(f"SMILES: {s}\nResult: {result}, Reason: {reason}\n")