"""
Classifies: CHEBI:48927 N-acyl-L-alpha-amino acid
"""
"""
Classifies: N-acyl-L-alpha-amino acid 
Definition: Any L-alpha-amino acid carrying an N-acyl substituent.
An L-alpha amino acid has a chiral alpha‐carbon bonded to an amino group and a carboxyl group.
The N-acyl substituent is defined by an acyl group (R–C(=O)–) attached to the amino nitrogen.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_N_acyl_L_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acyl L-alpha-amino acid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an N-acyl L-alpha-amino acid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for L-alpha-amino acid backbone.
    # We look for a chiral carbon with an amino group and a carboxyl group.
    aa_smarts1 = Chem.MolFromSmarts("[C@H](N)(C(=O)[O;H,-])")
    aa_smarts2 = Chem.MolFromSmarts("[C@@H](N)(C(=O)[O;H,-])")
    
    # Find matches for either chiral pattern.
    matches1 = mol.GetSubstructMatches(aa_smarts1)
    matches2 = mol.GetSubstructMatches(aa_smarts2)
    aa_matches = matches1 + matches2
    if not aa_matches:
        return False, "No L-alpha-amino acid backbone found"
    
    # Helper function to check whether a given carbon is part of an acyl group.
    def is_carbonyl_carbon(acyl_carbon):
        """
        Check if a carbon atom is part of a carbonyl group, i.e.
        it has exactly one double bond to an oxygen.
        """
        oxy_double = 0
        for bond in acyl_carbon.GetBonds():
            # Check that bond is a double bond and the neighbor is oxygen
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                nbr = bond.GetOtherAtom(acyl_carbon)
                if nbr.GetAtomicNum() == 8:
                    oxy_double += 1
        return oxy_double >= 1

    # Iterate over each amino acid backbone match and check for N-acyl substitution.
    # In the SMARTS pattern, the first atom (index 0) is the chiral alpha carbon.
    for match in aa_matches:
        alpha_idx = match[0]
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        
        # Get neighbors of the alpha carbon.
        neighbors = alpha_atom.GetNeighbors()
        
        # Identify the amino nitrogen: it should be a nitrogen atom.
        amino_nitrogens = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 7]
        if not amino_nitrogens:
            continue  # This candidate doesn't have an amino group; skip.
        
        # For each amino nitrogen, check if it is acylated.
        for amine in amino_nitrogens:
            # List substituents attached to the nitrogen that are not the alpha carbon.
            subs = [nbr for nbr in amine.GetNeighbors() if nbr.GetIdx() != alpha_idx]
            for sub in subs:
                # We are looking for a carbon that belongs to an acyl group.
                # Check if the substituent is carbon and has a carbonyl group feature.
                if sub.GetAtomicNum() == 6:
                    # Check if this carbon has a double bonded oxygen (i.e., is carbonyl).
                    if is_carbonyl_carbon(sub):
                        # Found an acyl group on the amino nitrogen.
                        return True, "Contains L-alpha-amino acid backbone with acylated amino group"
    
    return False, "Found L-alpha-amino acid backbone, but amino group is not acylated"
    
# If run as main, you might test the function with an example.
if __name__ == "__main__":
    # Example: N-acetyl-L-aspartic acid SMILES
    example_smiles = "CC(=O)N[C@@H](CC(O)=O)C(O)=O"
    result, reason = is_N_acyl_L_alpha_amino_acid(example_smiles)
    print("Result:", result)
    print("Reason:", reason)