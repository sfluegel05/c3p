"""
Classifies: CHEBI:86315 methyl sulfide
"""
"""
Classifies: Methyl sulfide
Definition: Any aliphatic sulfide in which at least one of the organyl groups attached to the sulfur is a methyl group.
An aliphatic sulfide (thioether) is an S atom connected only by single bonds to other atoms.
"""

from rdkit import Chem

def is_methyl_sulfide(smiles: str):
    """
    Determines if a molecule is a methyl sulfide based on its SMILES string.
    An aliphatic (thioether) sulfide must have an S atom only single-bonded to its substituents,
    and at least one of those substituents must be a methyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl sulfide, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate through all the atoms to find sulfur atoms (S)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:
            # Check that the S atom is only involved in single bonds;
            # if it has any non-single bond, it is not a typical thioether
            bonds = atom.GetBonds()
            if not all(bond.GetBondType() == Chem.rdchem.BondType.SINGLE for bond in bonds):
                continue  # Skip if S is not in a simple sulfide environment
            
            # Examine substituents (neighbors) of the sulfur atom.
            # We require at least one CH3 group, which will typically have atomic number 6 
            # and be bonded to only one heavy atom (the sulfur), because the hydrogens are implicit.
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetDegree() == 1:
                    return True, "Found a sulfur atom in a thioether-attached methyl group"
                    
    # If no S atom with a methyl substituent is found, return False.
    return False, "No aliphatic sulfide with a methyl substituent found"

# Example usage:
if __name__ == "__main__":
    test_smiles = "CSCC(O)=O"  # (methylthio)acetic acid example
    result, reason = is_methyl_sulfide(test_smiles)
    print(f"Result: {result}\nReason: {reason}")