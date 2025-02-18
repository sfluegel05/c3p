"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: Thiol – An organosulfur compound in which a thiol group (-SH) 
is attached to a carbon atom of any aliphatic or aromatic moiety.
However, when the -SH is part of an amino acid/peptide backbone (e.g. in cysteine
or peptide sequences), it should not be classified as a simple thiol.
"""

from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a (nonpeptidic) thiol based on its SMILES string.
    
    A thiol is defined as having at least one -SH group attached to a carbon where the 
    local environment does not suggest that the thiol is part of an amino acid or peptide 
    backbone. In particular, in many peptides (or free amino acids like cysteine/homocysteine)
    the -SH group is attached to a CH2 that is also (through its only neighbor) attached to a 
    chiral α-carbon that bears an amine and a carboxyl group.
    
    For each putative thiol group we:
      1. Look for a carbon attached to sulfur where the sulfur has exactly one hydrogen (when Hs are added).
      2. Check whether the CH2 (or CH) bonded to S is also attached to a carbon that looks like the 
         α–carbon of an amino acid (i.e. having a nitrogen and a carbonyl fragment connected).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if at least one -SH group is found on a carbon that is not in a peptide-like environment.
        str: Explanation / Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that the hydrogen on S becomes explicit.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS for a thiol group: a carbon (atomic number 6) bonded to an S that has exactly one H.
    # This is written as a bond from a carbon to an S with X2 (valence 2) that has one H.
    thiol_pattern = Chem.MolFromSmarts("[#6]-[S;X2&H1]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)
    if not thiol_matches:
        return False, "No thiol group (-SH attached to a carbon) found in the molecule"

    # For each matching thiol group, check if its environment seems peptide‐like.
    # Heuristic: for the carbon that is directly bonded to S (let's call it C_thiol),
    # check its neighbors (other than the sulfur) to see if any look like amino acid alpha carbons.
    # Typical alpha-carbon in amino acids is bonded to (at least) one nitrogen and one carbon that carries a C=O.
    for match in thiol_matches:
        carbon_idx, sulfur_idx = match  # match gives (index_of_C, index_of_S)
        C_thiol = mol.GetAtomWithIdx(carbon_idx)
        is_peptide_env = False   # flag for this thiol group
        
        # Examine neighbors of the thiol-bearing carbon (except the sulfur itself).
        for nb in C_thiol.GetNeighbors():
            if nb.GetIdx() == sulfur_idx:
                continue
            # Only consider if the neighbor is another carbon atom (which could be the alpha–carbon)
            if nb.GetAtomicNum() == 6:
                # Check if this neighbor (potential α–carbon) has at least one nitrogen neighbor ...
                has_nitrogen = any(nbr.GetAtomicNum() == 7 for nbr in nb.GetNeighbors())
                # ... and at least one neighboring carbon that is double-bonded to an oxygen.
                has_carbonyl = False
                for second_nb in nb.GetNeighbors():
                    # Skip if second_nb is C_thiol to avoid looping back.
                    if second_nb.GetIdx() == C_thiol.GetIdx():
                        continue
                    if second_nb.GetAtomicNum() == 6:
                        # Check bonds on this carbon for a double bond to oxygen.
                        for bond in second_nb.GetBonds():
                            # bond order 2.0 indicates a double bond.
                            if bond.GetBondTypeAsDouble() == 2.0:
                                other = bond.GetOtherAtom(second_nb)
                                if other.GetAtomicNum() == 8:
                                    has_carbonyl = True
                                    break
                    if has_carbonyl:
                        break
                if has_nitrogen and has_carbonyl:
                    is_peptide_env = True
                    break  # one peptide-like environment is enough
                
        # If we found at least one thiol that is NOT in a peptide environment, classify as thiol.
        if not is_peptide_env:
            return True, "Molecule contains a thiol group (-SH) attached to a carbon atom outside a peptide backbone"
    
    # If every detected thiol group seems to be tied to an amino acid/peptide-like environment, do not classify.
    return False, "Thiols detected are associated with peptide backbones or amino acid environments, not classified as simple thiol compounds"


# Example usage (for testing purposes)
if __name__ == "__main__":
    # A few examples:
    test_molecules = {
        "2-Methoxybenzenethiol": "SC=1C(OC)=CC=CC1",
        "Cysteamine hydrochloride": "Cl.SCCN",
        "3-mercaptopyruvic acid": "OC(=O)C(=O)CS",
        "Arg-Cys-Asp (peptide example)": "SC[C@H](NC(=O)[C@@H](N)CCCN=C(N)N)C(=O)N[C@@H](CC(O)=O)C(O)=O"
    }
    
    for name, smi in test_molecules.items():
        result, reason = is_thiol(smi)
        print(f"{name}:\n  SMILES: {smi}\n  Classified as thiol? {result}\n  Reason: {reason}\n")