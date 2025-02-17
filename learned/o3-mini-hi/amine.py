"""
Classifies: CHEBI:32952 amine
"""
"""
Classifies: Amine
Definition: A compound formally derived from ammonia by replacing one, two or three hydrogen atoms
by hydrocarbyl groups.
"""

from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    A compound qualifies as an amine if it contains at least one nitrogen atom that is
    derived from ammonia by replacing one, two or three hydrogens with hydrocarbyl groups.
    (In other words, the nitrogen should have three substituents when counting implicit hydrogens,
    and should not be part of an amide function.)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as having an amine functional group; False otherwise
        str: Reason for the classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS for an amine:
    # [#7] represents any nitrogen atom.
    # X3 makes sure it has three neighbors when implicit hydrogens are counted (thus primary, secondary,
    # or tertiary amine).
    # The substructure "!$(#7-C(=O))" excludes cases where the nitrogen is directly bound to a carbonyl
    # (e.g. an amide), since amides are not considered amines in this context.
    amine_smarts = "[#7;X3;!$(#7-C(=O))]"
    amine_pattern = Chem.MolFromSmarts(amine_smarts)
    if amine_pattern is None:
        return False, "Error in the SMARTS pattern for detecting amines"
    
    # Find all substructure matches of the defined amine pattern within the molecule.
    matches = mol.GetSubstructMatches(amine_pattern)
    
    if matches:
        return True, f"Found {len(matches)} amine group(s) that are not part of amide functionalities."
    else:
        return False, "No amine functional group found (or present nitrogen atoms are in non-amine environments)."