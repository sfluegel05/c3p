"""
Classifies: CHEBI:87659 dodecanoate ester
"""
"""
Classifies: CHEBI:36000 dodecanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dodecanoate_ester(smiles: str):
    """
    Determines if a molecule is a dodecanoate ester based on its SMILES string.
    A dodecanoate ester is an ester where the acyl chain (acid component) is derived from lauric acid 
    (dodecanoic acid, a saturated 12-carbon fatty acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dodecanoate ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define lauric acid (dodecanoic acid)
    lauric_acid_smiles = 'CCCCCCCCCCC(=O)O'  # 12 carbons including carbonyl carbon
    lauric_acid_mol = Chem.MolFromSmiles(lauric_acid_smiles)
    lauric_acid_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(lauric_acid_mol, radius=2)

    # Identify ester groups: pattern C(=O)O
    ester_pattern = Chem.MolFromSmarts('C(=O)O')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester groups found"

    for match in ester_matches:
        carbonyl_c_idx = match[0]  # Carbonyl carbon
        carbonyl_o_idx = match[1]  # Carbonyl oxygen
        ester_o_idx = match[2]     # Ester oxygen

        # Check bond between carbonyl carbon and ester oxygen
        bond = mol.GetBondBetweenAtoms(carbonyl_c_idx, ester_o_idx)
        if bond is None:
            continue  # Skip if bond does not exist

        bond_to_break = bond.GetIdx()

        # Break the bond between carbonyl carbon and ester oxygen
        fragmented_mol = Chem.FragmentOnBonds(mol, [bond_to_break], addDummies=False)

        # Get the fragments
        frags = Chem.GetMolFrags(fragmented_mol, asMols=True, sanitizeFrags=True)
        if not frags:
            continue  # No fragments obtained, skip to next match

        # The acyl fragment should contain the carbonyl carbon
        acyl_frag = None
        for frag in frags:
            atom_indices = [atom.GetIdx() for atom in frag.GetAtoms()]
            if carbonyl_c_idx in atom_indices:
                acyl_frag = frag
                break
        if acyl_frag is None:
            continue  # Could not find acyl fragment, continue to next ester group

        # Sanitize the fragment to fix valencies
        Chem.SanitizeMol(acyl_frag)

        # Compare the acyl fragment to lauric acid
        # Compute fingerprints
        acyl_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(acyl_frag, radius=2)
        # Calculate Tanimoto similarity
        similarity = Chem.DataStructs.TanimotoSimilarity(acyl_fp, lauric_acid_fp)
        if similarity > 0.9:
            return True, "Contains dodecanoate ester group"

        # Alternatively, compare molecular formulas
        lauric_formula = rdMolDescriptors.CalcMolFormula(lauric_acid_mol)
        acyl_formula = rdMolDescriptors.CalcMolFormula(acyl_frag)
        if lauric_formula == acyl_formula:
            return True, "Contains dodecanoate ester group based on molecular formula"

    return False, "No dodecanoate ester groups with lauric acid found"