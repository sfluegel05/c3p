"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: Thiol – An organosulfur compound in which a thiol group (-SH) 
is attached to a carbon atom of any aliphatic or aromatic moiety.
Exceptions: When the S-containing group is part of an amino acid or peptide backbone 
(i.e. the carbon attached to S appears as an α–carbon with both an amine and a carbonyl-bearing neighbor),
the molecule should not be classified as a simple thiol.
"""

from rdkit import Chem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a nonpeptidic thiol based on its SMILES string.
    
    Strategy:
      1. Parse the molecule and add implicit hydrogens.
      2. Loop over all sulfur atoms (atomic number 16). For each S:
           a. Skip it if it is oxidized (i.e. S is double-bonded to an oxygen).
           b. Check that the S atom has at least one hydrogen attached (for an –SH group).
           c. Check that the S atom is attached to exactly one heavy (non-hydrogen) neighbor.
              (A genuine thiol –SH should have one attachment to a carbon.)
           d. Ensure that the attached carbon is not in a peptide-like environment.
              We use a heuristic: if this carbon (C_thiol) has at least one nitrogen neighbor and 
              at least one neighbor (other than the S) that is a carbon double-bonded to oxygen (a carbonyl), 
              then we consider it peptide-like.
      3. If any sulfur passing these tests is found, return True with an explanation.
      4. Otherwise return False.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a valid (nonpeptidic) thiol group is found, False otherwise.
        str: Explanation / Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For our purposes the implicit hydrogens are used via GetTotalNumHs().
    # Loop over every atom to find S (atomic number 16)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue  # not sulfur
        
        # Skip oxidized sulfur atoms (i.e. S with a double bond to O).
        if any(bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(atom).GetAtomicNum() == 8
               for bond in atom.GetBonds()):
            continue
        
        # Check that the sulfur has at least one hydrogen attached.
        # In a genuine thiol group (-SH), the S should carry one hydrogen.
        if atom.GetTotalNumHs() < 1:
            continue
        
        # Get heavy (non-hydrogen) neighbors
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # For a free thiol (-SH) the S should be attached to exactly one heavy atom.
        if len(heavy_neighbors) != 1:
            continue
        
        neighbor = heavy_neighbors[0]
        # The heavy neighbor must be a carbon to qualify as a typical thiol.
        if neighbor.GetAtomicNum() != 6:
            continue
        
        # Heuristic to detect peptide-like environment:
        # Check if the carbon (C_thiol) attached to S has:
        #   a) at least one neighboring nitrogen, and 
        #   b) at least one neighboring carbon (other than the S) that is double-bonded to an oxygen.
        peptide_like = False
        neighbors_of_c = [n for n in neighbor.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
        
        has_nitrogen = any(n.GetAtomicNum() == 7 for n in neighbors_of_c)
        
        has_carbonyl = False
        for n in neighbors_of_c:
            if n.GetAtomicNum() == 6:
                # Look for a double bond from this carbon to an oxygen.
                for bond in n.GetBonds():
                    # bond.GetOtherAtom(n) returns the neighbor
                    if bond.GetBondTypeAsDouble() == 2.0 and bond.GetOtherAtom(n).GetAtomicNum() == 8:
                        has_carbonyl = True
                        break
            if has_carbonyl:
                break
        
        if has_nitrogen and has_carbonyl:
            peptide_like = True
        
        # If this S--C link is not peptide-like, we count it as a valid thiol group.
        if not peptide_like:
            return True, "Molecule contains a free thiol (-SH) group attached to a carbon outside a peptide backbone."
    
    return False, "No valid (nonpeptidic) thiol (-SH) group found."

# Example usage (for testing purposes)
if __name__ == "__main__":
    test_molecules = {
        "2-methoxy-4-(methylthio)-N-pyridin-4-ylbenzamide": "COC1=C(C=CC(=C1)SC)C(=O)NC2=CC=NC=C2",
        "(S)-fluoxytioconazole (thiol)": "O[C@](CN1N=CN=C1S)(C1=CC=C(F)C=C1F)C(F)(F)C1=CC=C(OC2=CC=C(C=C2)C#N)C=N1",
        "coenzyme A": "CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12",
        "3-mercaptopyruvic acid": "OC(=O)C(=O)CS",
        "Cysteamine hydrochloride": "Cl.SCCN",
        "mercaptoethanol": "OCCS",
        "cysteamine": "NCCS",  # should be classified as thiol
    }
    
    for name, smi in test_molecules.items():
        result, reason = is_thiol(smi)
        print(f"{name}:\n  SMILES: {smi}\n  Classified as thiol? {result}\n  Reason: {reason}\n")