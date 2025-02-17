"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: lactol 
Definition: Cyclic hemiacetals formed by intramolecular addition of a hydroxy group 
to an aldehydic or ketonic carbonyl group. In a lactol the hemiacetal (anomeric) carbon is sp3,
not a carbonyl, and is bonded to one –OH group and one –OR group (where the –OR belongs to a ring).
"""

from rdkit import Chem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    
    Heuristic approach:
      1. Parse the SMILES and add explicit hydrogens.
      2. Define a SMARTS pattern that should match a cyclic hemiacetal center. 
         This is the pattern: [C;R;!$(C=O)]([OX2H])([O;R])
         • [C;R;!$(C=O)]: a carbon atom (in a ring) that is not double‐bonded to oxygen
         • ([OX2H]): one substituent is an –OH group (an oxygen with at least one hydrogen)
         • ([O;R]): the other substituent is an oxygen that belongs to a ring (i.e. an ether oxygen in a ring)
      3. If any match is found the molecule is classified as a lactol.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if a lactol center is identified, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES to an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Add explicit hydrogens to ensure –OH groups are recognized
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a cyclic hemiacetal (lactol) center.
    # This matches a ring carbon that is not a carbonyl and has exactly one attached –OH and one attached ether oxygen
    lactol_smarts = "[C;R;!$(C=O)]([OX2H])([O;R])"
    lactol_pattern = Chem.MolFromSmarts(lactol_smarts)
    if lactol_pattern is None:
        return False, "Error creating SMARTS pattern"
    
    # Search for the lactol pattern in the molecule
    matches = mol.GetSubstructMatches(lactol_pattern)
    if matches:
        # If one or more matches are found, we assume the molecule contains a lactol center.
        # Report the first matching atom index as the lactol center.
        lactol_carbon_idx = matches[0][0]
        return True, f"Found cyclic hemiacetal center at carbon atom index {lactol_carbon_idx}"
    
    return False, "No cyclic hemiacetal (lactol) center found"

# Example usage (optional):
if __name__ == "__main__":
    # Test example: beta-ascarylopyranose is known to be a lactol.
    smiles_example = "C[C@@H]1O[C@H](O)[C@H](O)C[C@H]1O"
    result, reason = is_lactol(smiles_example)
    print(result, reason)