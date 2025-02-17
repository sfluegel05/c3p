"""
Classifies: CHEBI:140310 phenyl acetates
"""
"""
Classifies: Phenyl Acetates
Definition: An acetate ester obtained by formal condensation of the carboxy group of acetic acid with 
the hydroxy group of any phenol.
"""

from rdkit import Chem

def is_phenyl_acetates(smiles: str):
    """
    Determines if a molecule is a phenyl acetate based on its SMILES string.
    
    A phenyl acetate has an aromatic ring (typically a phenol) which has been acetylated,
    forming an ester bond: -O-C(=O)CH3. We look for such a pattern where the oxygen is directly 
    attached to an aromatic (c) atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains a phenyl acetate moiety, False otherwise.
        str: Explanation for the classification.
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for the phenyl acetate group.
    # This pattern looks for an aromatic carbon (c) directly bonded to an oxygen (O)
    # which is further bonded to an acyl group C(=O)C (i.e. acetyl group).
    acetate_pattern = Chem.MolFromSmarts("cOC(=O)C")
    if acetate_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Check if the pattern is part of the molecule.
    if mol.HasSubstructMatch(acetate_pattern):
        return True, "Molecule contains a phenyl acetate moiety (acetylated phenolic oxygen)."
    else:
        return False, "No phenyl acetate substructure (cOC(=O)C) found in the molecule."

# Example usage:
if __name__ == "__main__":
    test_smiles = "CC(=O)Oc1ccccc1"  # phenyl acetate
    result, reason = is_phenyl_acetates(test_smiles)
    print(result, reason)