"""
Classifies: CHEBI:23053 catechin
"""
"""
Classifies: CHEBI:23053 catechin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechin(smiles: str):
    """
    Determines if a molecule is a catechin based on its SMILES string.
    A catechin is a member of the class of hydroxyflavan that has a flavan-3-ol skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define core structure pattern (flavan-3-ol skeleton)
    # This pattern matches the basic 2-phenylchromane structure with at least one hydroxyl group
    # and allows for various substitutions
    core_pattern = Chem.MolFromSmarts("C1C(O)Cc2c(O)cc(O)cc2O1")
    if not mol.HasSubstructMatch(core_pattern):
        return False, "No flavan-3-ol core structure found"

    # Check for characteristic features
    # 1. At least two aromatic rings (one from flavan, one from phenyl)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    if aromatic_rings < 2:
        return False, f"Found {aromatic_rings} aromatic rings, need at least 2"

    # 2. Check for hydroxyl groups (both free and in substitutions)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = len(mol.GetSubstructMatches(hydroxyl_pattern))
    
    # 3. Check for common substitutions
    substitution_patterns = [
        Chem.MolFromSmarts("[OH]c1cc(O)c(O)c(O)c1"),  # Gallate
        Chem.MolFromSmarts("[OH]S(=O)(=O)[OH]"),      # Sulfate
        Chem.MolFromSmarts("C(=O)O"),                 # Ester
        Chem.MolFromSmarts("C(=O)[OH]"),              # Carboxyl
        Chem.MolFromSmarts("C=C"),                    # Double bond
        Chem.MolFromSmarts("C#N"),                    # Nitrile
        Chem.MolFromSmarts("C1OC1"),                  # Epoxide
        Chem.MolFromSmarts("CO"),                     # Methoxy
        Chem.MolFromSmarts("C=O"),                    # Carbonyl
    ]
    
    has_substitution = any(mol.HasSubstructMatch(patt) for patt in substitution_patterns)
    
    # Adjust hydroxyl count requirements based on substitutions
    if has_substitution:
        min_hydroxyl = 2  # Allow fewer hydroxyls if substitutions are present
    else:
        min_hydroxyl = 3  # Basic catechin needs at least 3 hydroxyls

    if hydroxyl_matches < min_hydroxyl:
        return False, f"Found {hydroxyl_matches} hydroxyl groups, need at least {min_hydroxyl}"

    # 4. Check molecular weight range (200-1000 Da typical for catechin derivatives)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 1000:
        return False, f"Molecular weight {mol_wt:.1f} outside typical catechin range"

    # 5. Additional checks for common catechin features
    # Check for at least one aromatic hydroxyl group
    aromatic_oh_pattern = Chem.MolFromSmarts("c[OH]")
    if not mol.HasSubstructMatch(aromatic_oh_pattern):
        return False, "No aromatic hydroxyl groups found"

    return True, "Contains flavan-3-ol core structure with characteristic features"