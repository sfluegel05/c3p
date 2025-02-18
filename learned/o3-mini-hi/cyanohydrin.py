"""
Classifies: CHEBI:23437 cyanohydrin
"""
"""
Classifies: Cyanohydrin (alpha-hydroxynitrile)
An alpha-hydroxynitrile results from the formal addition of hydrogen cyanide 
to the carbonyl (C=O) bond of an aldehyde or ketone, yielding a molecule with a 
hydroxyl and a nitrile group on the same carbon.
"""
from rdkit import Chem

def is_cyanohydrin(smiles: str):
    """
    Determines if a molecule is a cyanohydrin based on its SMILES string.
    A cyanohydrin contains a sp3-hybridized carbon with both an -OH group and a nitrile (-C≡N) attached.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a cyanohydrin, False otherwise.
        str: A reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Compile SMARTS pattern for a cyanohydrin group:
    # Explanation:
    #   [#6X4]: a sp3 (tetrahedral) carbon,
    #   ([OX2H]): attached to an oxygen with H (hydroxyl),
    #   ([#6X2]#[#7]): attached to a sp-hybridized carbon triple bonded to a nitrogen (nitrile).
    cyanohydrin_smarts = "[#6X4]([OX2H])([#6X2]#[#7])"
    pattern = Chem.MolFromSmarts(cyanohydrin_smarts)
    
    if pattern is None:
        # This should not happen; SMARTS string is expected to be valid.
        return False, "Invalid SMARTS pattern for cyanohydrin"
    
    # Check if the molecule has at least one match for the cyanohydrin pattern.
    if mol.HasSubstructMatch(pattern):
        return True, "Cyanohydrin pattern (carbon with both -OH and -C≡N substituents) detected"
    else:
        return False, "Cyanohydrin pattern not found in the molecule"