"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: Thiol – An organosulfur compound in which a thiol group (-SH) is attached 
to a carbon atom of any aliphatic or aromatic moiety. However, when the S-containing group 
is part of a peptide backbone (i.e. the S is attached to a peptide α–carbon which has both 
an amine and a carbonyl neighbor), the molecule should not be classified as a simple thiol.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_thiol(smiles: str):
    """
    Determines if a molecule is a valid (nonpeptidic) thiol based on its SMILES string.
    
    Strategy:
      1. Parse the molecule.
      2. Loop over all sulfur atoms (atomic number 16). For each S:
           a. Skip if the sulfur is oxidized (i.e. attached via a double bond to oxygen).
           b. Check that the S atom has at least one hydrogen attached (for an –SH group).
           c. Require that S is attached to exactly one heavy (non-hydrogen) neighbor.
           d. Check that the heavy neighbor is a carbon.
           e. If that carbon is aromatic, accept the group.
           f. Otherwise (for an aliphatic carbon) we check if that carbon looks peptide‐like:
              If the carbon (excluding the S neighbor) has at least one nitrogen neighbor 
              and at least one neighbor that is a carbon double-bonded to an oxygen (a carbonyl), 
              then we mark it as peptide-like and skip this S.
      3. If any S passes these tests, return True with an explanation.
         Otherwise return False with the reason.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule contains a nonpeptidic free thiol (-SH) group, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Loop over every atom to check for sulfur atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 16:
            continue  # not a sulfur
        
        # Skip S if it is oxidized (i.e. S has a double bond to an oxygen)
        oxidized = False
        for bond in atom.GetBonds():
            # Use GetBondType() for a more robust double-bond check
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                other = bond.GetOtherAtom(atom)
                if other.GetAtomicNum() == 8:
                    oxidized = True
                    break
        if oxidized:
            continue
        
        # Check that the sulfur atom has at least one hydrogen attached (free –SH group) 
        # (implicit hydrogens are taken into account)
        if atom.GetTotalNumHs() < 1:
            continue
        
        # Get the heavy (non-hydrogen) neighbors of S
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # A free thiol (-SH) should be attached to exactly one heavy atom.
        if len(heavy_neighbors) != 1:
            continue
        
        neighbor = heavy_neighbors[0]
        # The neighbor must be a carbon (atomic number 6) for thiol classification.
        if neighbor.GetAtomicNum() != 6:
            continue

        # If the attached carbon is aromatic, we assume it is a valid thiol group.
        if neighbor.GetIsAromatic():
            return True, "Molecule contains a free thiol (-SH) group attached to an aromatic carbon."
        
        # For an aliphatic carbon, check if it appears peptide-like.
        # Heuristic: if the carbon (excluding the S) has at least one nitrogen neighbor
        # and at least one neighbor (other than S) which is a carbon that is double-bonded to O,
        # then assume it is part of a peptide backbone.
        neighbors_of_c = [n for n in neighbor.GetNeighbors() if n.GetIdx() != atom.GetIdx()]
        has_nitrogen = any(n.GetAtomicNum() == 7 for n in neighbors_of_c)
        has_carbonyl = False
        for n in neighbors_of_c:
            if n.GetAtomicNum() == 6:
                for bond in n.GetBonds():
                    if bond.GetBondType() == rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(n)
                        if other.GetAtomicNum() == 8:
                            has_carbonyl = True
                            break
            if has_carbonyl:
                break
        
        if has_nitrogen and has_carbonyl:
            # The S–C link appears to be in a peptide-like environment; skip it.
            continue
        else:
            return True, ("Molecule contains a free thiol (-SH) group attached to a carbon "
                          "outside a peptide backbone.")
    
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