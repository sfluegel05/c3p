"""
Classifies: CHEBI:17334 penicillin
"""
#!/usr/bin/env python
"""
Classifies: Penicillins
Definition:
  “Any member of the group of substituted penams containing two methyl substituents at position 2,
   a carboxylate substituent at position 3 and a carboxamido group at position 6.”
   
This implementation uses heuristic SMARTS patterns:
  - A relaxed SMARTS pattern for a fused bicyclic penam core (β‐lactam fused to thiazolidine)
  - A pattern to detect two methyl substituents attached to the S-bound carbon (position 2)
  - A generic carboxylate fragment (C(=O)O or C(=O)[O-])
  - A generic carboxamido fragment (N-C(=O))
Note:
  Penicillin structures are variable and may contain extra stereochemistry or substituents.
  We ignore stereochemistry in the core SMARTS.
"""

from rdkit import Chem

def is_penicillin(smiles: str):
    """
    Determines if a molecule (given as a SMILES string) is a penicillin.
    A molecule is considered a penicillin if it satisfies these conditions:
      - Contains a fused bicyclic penam core (β‐lactam fused to thiazolidine)
      - Contains two methyl substituents on the carbon adjacent to sulfur (position 2)
      - Contains a carboxylate substituent (e.g. C(=O)O or C(=O)[O-])
      - Contains a carboxamido group (N-C(=O))
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a penicillin, otherwise False.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Strategy improvement:
    # We relax the penam core SMARTS by removing stereochemical descriptors.
    # The pattern below tries to capture a fused bicyclic system having:
    #   - A four-membered beta-lactam (with an amide carbonyl)
    #   - A five-membered thiazolidine containing a sulfur and a carbon that bears two methyl groups.
    #
    # The pattern below does not include chirality and is intended to be a looser match.
    penam_core_smarts = "N1C(=O)C2SC(C)(C)CN2C1=O"
    penam_core = Chem.MolFromSmarts(penam_core_smarts)
    if penam_core is None:
        return False, "Error in penam core SMARTS"
    
    if not mol.HasSubstructMatch(penam_core):
        return False, "Molecule does not contain the expected penam core (relaxed match)"
    
    # Check for two methyl substituents at the key carbon.
    # We use a simple SMARTS that matches an S-bound carbon that has two CH3 groups.
    # This pattern is less strict regarding bond direction and stereochemistry.
    dimethyl_pattern = Chem.MolFromSmarts("S-C([CH3])([CH3])")
    if dimethyl_pattern is None:
        return False, "Error in dimethyl SMARTS"
    if not mol.HasSubstructMatch(dimethyl_pattern):
        return False, "Missing two methyl substituents at position 2"
    
    # Check for a carboxylate substituent.
    # We match either the deprotonated form (C(=O)[O-]) or the acid form (C(=O)O)
    carboxylate1 = Chem.MolFromSmarts("C(=O)[O-]")
    carboxylate2 = Chem.MolFromSmarts("C(=O)O")
    if carboxylate1 is None or carboxylate2 is None:
        return False, "Error in carboxylate SMARTS"
    if not (mol.HasSubstructMatch(carboxylate1) or mol.HasSubstructMatch(carboxylate2)):
        return False, "Missing a carboxylate substituent (C(=O)O group)"
    
    # Check for a carboxamido group by looking for an N-C(=O) fragment.
    carboxamido = Chem.MolFromSmarts("N-C(=O)")
    if carboxamido is None:
        return False, "Error in carboxamido SMARTS"
    if not mol.HasSubstructMatch(carboxamido):
        return False, "Missing a carboxamido group (N-C(=O) fragment)"
    
    return True, "Molecule has a penam core with dimethyl at position 2, a carboxylate at position 3, and a carboxamido group at position 6."

# Example usage if run as a script:
if __name__ == "__main__":
    # Test with one of the provided penicillin SMILES (e.g., Penicillin K)
    test_smiles = "CCCCCCCC(=O)N[C@H]1[C@H]2SC(C)(C)[C@@H](N2C1=O)C(O)=O"
    is_pen, reason = is_penicillin(test_smiles)
    print("Is penicillin?", is_pen)
    print("Reason:", reason)