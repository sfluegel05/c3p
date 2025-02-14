"""
Classifies: CHEBI:50477 butyrate ester
"""
"""
Classifies: butyrate ester

Definition: Any carboxylic ester where the carboxylic acid component is butyric acid.
"""

from rdkit import Chem

def is_butyrate_ester(smiles: str):
    """
    Determines if a molecule is a butyrate ester based on its SMILES string.
    A butyrate ester is any carboxylic ester where the carboxylic acid component is butyric acid.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a butyrate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for ester group
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;!$(C(=O))]")
    if ester_pattern is None:
        return False, "Failed to create ester pattern"

    # Find ester groups
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    # Define butyric acid molecule (butyric acid)
    butyric_acid = Chem.MolFromSmiles("CCCC(=O)O")
    butyric_acid_pattern = Chem.MolFromSmarts("CCCC(=O)[O;H]")

    # Define butyryl group (butyric acid minus the OH)
    butyryl_group = Chem.MolFromSmiles("CCCC(=O)")
    if butyryl_group is None:
        return False, "Failed to create butyryl group pattern"

    # Loop over ester groups
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon index
        ester_o_idx = match[1]     # Ester oxygen index
        alcohol_c_idx = match[2]   # Alcohol carbon index

        # Break the bond between ester oxygen and alcohol carbon to get acyl fragment
        bond = mol.GetBondBetweenAtoms(ester_o_idx, alcohol_c_idx)
        if bond is None:
            continue
        bond_idx = bond.GetIdx()

        # Fragment the molecule at the ester bond
        fragments = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
        # Get the acyl fragment (should contain the carbonyl carbon)
        frags = Chem.GetMolFrags(fragments, asMols=True, sanitizeFrags=True)
        acyl_frag = None
        for frag in frags:
            # Look for the carbonyl carbon connected to [*][C](=O)[O][*]
            if frag.HasSubstructMatch(Chem.MolFromSmarts("[#6](=O)[O][#1]")):
                acyl_frag = frag
                break
        if acyl_frag is None:
            continue

        # Check if acyl fragment matches butyric acid minus the OH group
        if acyl_frag.HasSubstructMatch(butyryl_group):
            return True, "Contains butyrate ester group"

    return False, "Does not contain butyrate ester group"