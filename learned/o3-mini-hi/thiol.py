"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: Thiol – An organosulfur compound in which a thiol group (-SH) is attached 
to a carbon atom (aliphatic or aromatic). 
Note: To avoid mis‐classifying peptides (which may contain –SH groups on cysteine residues)
as simple thiol compounds, we apply an additional check to see if the –SH group is part of a 
peptide-like backbone.
"""
from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a thiol (nonpeptidic) based on its SMILES string.
    A thiol is defined as an organosulfur compound containing at least one -SH group attached 
    to a carbon atom that is not clearly part of a peptide backbone.
    
    The strategy is first to add explicit hydrogens and then search for a SMARTS pattern for a 
    thiol group defined as a carbon atom bonded to a sulfur atom that has exactly one hydrogen.
    Then, for each such thiol group found, we check whether the carbon (to which the sulfur is attached)
    is also attached to a nitrogen that in turn is bonded to a carbonyl group – a pattern characteristic 
    of the peptide backbone. If all observed thiol groups are part of such a peptide pattern, 
    the molecule is not classified as a (simple) thiol.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule has at least one non‐peptide thiol group, else False.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
 
    # Add hydrogens so that the hydrogen on sulfur is explicitly represented.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a thiol group: a carbon (#6) bonded to a sulfur atom that
    # has exactly one hydrogen and a total valence of 2 (i.e. one bond to C and one to H).
    # This pattern should capture -SH groups attached to any carbon.
    thiol_pattern = Chem.MolFromSmarts("[#6]-[S;X2&H1]")
    
    # Look for matches of the thiol pattern in the molecule.
    matches = mol.GetSubstructMatches(thiol_pattern)
    if not matches:
        return False, "No thiol group (-SH attached to a carbon) found in the molecule"
        
    # For each thiol match, check if the carbon directly bonded to S looks like it is part of a peptide.
    # The heuristic: Does the carbon (attached to the sulfur) have a neighbor that is nitrogen
    # which in turn is bound to a carbonyl carbon? (i.e. a common peptide-bond motif: -C-CH2-S vs.
    # -C(N)C(=O)- which is typical in amino acids.)
    for match in matches:
        carbon_idx, sulfur_idx = match  # match returns a tuple (idx of C, idx of S)
        c_atom = mol.GetAtomWithIdx(carbon_idx)
        is_peptide_thiol = False
        for neighbor in c_atom.GetNeighbors():
            # We skip the atom that is the sulfur of the thiol.
            if neighbor.GetIdx() == sulfur_idx:
                continue
            # Look for a nitrogen neighbor.
            if neighbor.GetAtomicNum() == 7:
                # For each nitrogen bonded to our carbon, see if it is attached to a carbonyl.
                for nn in neighbor.GetNeighbors():
                    # Skip if neighbor is already our carbon.
                    if nn.GetIdx() == c_atom.GetIdx():
                        continue
                    # Check if this neighbor is carbon (atomic num 6)
                    if nn.GetAtomicNum() == 6:
                        # Examine bonds of this carbon to see if one is a double bond to oxygen.
                        for bond in nn.GetBonds():
                            # bond.GetBondTypeAsDouble() returns 2.0 for a double bond.
                            if bond.GetBondTypeAsDouble() == 2.0:
                                # Check if the other atom in the bond is oxygen.
                                other = bond.GetOtherAtom(nn)
                                if other.GetAtomicNum() == 8:
                                    # We found a N--C(=O) fragment; mark as peptide thiol.
                                    is_peptide_thiol = True
                                    break
                        if is_peptide_thiol:
                            break
                if is_peptide_thiol:
                    break
        # If we found at least one thiol group that is not part of a peptide backbone, we classify as thiol.
        if not is_peptide_thiol:
            return True, "Molecule contains a thiol group (-SH) attached to a carbon atom outside a peptide backbone"
    
    # If we arrive here, every -SH group found appears to be embedded in a peptide-like environment.
    return False, "Thiols detected are associated with peptide backbones, not classified as simple thiol compounds"
    
# Example usage (for testing purposes)
if __name__ == "__main__":
    # Testing with an example molecule: 2-Methoxybenzenethiol which should be classified as thiol.
    test_smiles = "SC=1C(OC)=CC=CC1"
    result, reason = is_thiol(test_smiles)
    print("Result:", result)
    print("Reason:", reason)