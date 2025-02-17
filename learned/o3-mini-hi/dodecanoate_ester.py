"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: any fatty acid ester in which the carboxylic acid component is lauric acid,
i.e. a dodecanoate ester. The acyl group must be exactly CH3(CH2)10C(=O)O.
This program builds a SMARTS that encodes:
    [CH3;H3]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-C(=O)O
so that if a match is found the laurate (dodecanoate) ester substructure exists,
and extra bonds on the “terminal” atoms would lower the hydrogen counts so that the SMARTS does not match.
"""

from rdkit import Chem

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule contains a dodecanoate (laurate) ester group.
    A dodecanoate ester here is defined as having the acyl moiety CH3(CH2)10C(=O)O.
    
    The procedure is as follows:
      1. Convert the input SMILES into an RDKit molecule.
      2. Use a SMARTS pattern that demands a terminal CH3 (with three implicit/explcit H’s)
         followed by exactly ten CH2 groups (with exactly two hydrogens each), finishing in a carbonyl group bound to an -O.
         This pattern does not match if the chain is extended or branched.
      3. If at least one substructure match is found then we classify the molecule as containing
         a laurate ester group.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule contains a laurate ester group, False otherwise.
        str: Explanation for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Build a SMARTS for the exact laurate ester acyl fragment:
    # CH3(CH2)10C(=O)O
    # Note: each atom is annotated with a hydrogen count constraint:
    #  - The terminal CH3 must have exactly 3 hydrogen atoms.
    #  - Each CH2 must have exactly 2 hydrogens.
    #  - Then the carbonyl carbon (without H count constraint) is bound to =O and O.
    # The SMARTS string written out is:
    #   [CH3;H3]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-[CH2;H2]-C(=O)O
    smarts = ("[CH3;H3]"
              "-[CH2;H2]"
              "-[CH2;H2]"
              "-[CH2;H2]"
              "-[CH2;H2]"
              "-[CH2;H2]"
              "-[CH2;H2]"
              "-[CH2;H2]"
              "-[CH2;H2]"
              "-[CH2;H2]"
              "-[CH2;H2]"
              "-C(=O)O")
              
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        return False, "Could not build SMARTS pattern for laurate ester."
    
    # Find matches of the laurate ester pattern in the molecule.
    matches = mol.GetSubstructMatches(query)
    if matches:
        return True, "Molecule contains a dodecanoate (laurate) ester group."
    else:
        # No match found – either the ester group is absent or 
        # the acyl chain is not exactly laurate (e.g. if it is extended or branched).
        return False, "No laurate ester group with an exact CH3(CH2)10C(=O)O acyl chain found."

# Example usage:
if __name__ == "__main__":
    # Test with a known laurate ester: 1-lauroyl-sn-glycerol.
    test_smiles = "CCCCCCCCCCCC(=O)OC[C@@H](O)CO"
    result, reasoning = is_dodecanoate_ester(test_smiles)
    print(result, reasoning)