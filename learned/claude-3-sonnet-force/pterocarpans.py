"""
Classifies: CHEBI:26377 pterocarpans
"""
"""
Classifies: CHEBI:51811 pterocarpan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_pterocarpans(smiles: str):
    """
    Determines if a molecule is a pterocarpan based on its SMILES string.
    Pterocarpans are members of the class of benzofurochromene with a
    6a,11a-dihydro-6H-[1]benzofuro[3,2-c]chromene skeleton and its substituted derivatives.
    They generally bear structural resemblance to isoflavanoids and are produced by plant tissues
    in response to infection. They are the 3,4-dihydroderivatives of coumestans.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pterocarpan, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pterocarpan scaffold pattern
    scaffold_pattern = Chem.MolFromSmarts("[C@]12Oc3ccccc3[C@@]1(O)COc1ccccc21")
    if not mol.HasSubstructMatch(scaffold_pattern):
        return False, "Does not contain the pterocarpan scaffold"

    # Check for fused rings
    fused_rings = AllChem.GetSymmSSSR(mol)
    if len(fused_rings) < 4:
        return False, "Does not have enough fused rings"

    # Check for required heteroatoms
    atom_counts = mol.GetAtomSMARTSPatternCount("O", 2, 4)
    if atom_counts <= 0:
        return False, "Does not contain the required oxygen atoms"

    # Check for common substitutions
    subs_pattern = Chem.MolFromSmarts("[OC,OCC]")
    subs_matches = mol.GetSubstructMatches(subs_pattern)
    if len(subs_matches) < 2:
        return False, "Does not have enough common substitutions"

    # Check for connectivity and size
    if not mol.GetRingInfo().AtomRings():
        return False, "Does not contain rings"
    if mol.GetNumHeavyAtoms() < 15:
        return False, "Too small to be a pterocarpan"

    # Check for structural similarity to isoflavanoids
    isoflavanoid_pattern = Chem.MolFromSmarts("[C@]12Oc3ccccc3[C@@]1(O)COc1ccccc21")
    if not mol.HasSubstructMatch(isoflavanoid_pattern):
        return False, "Does not resemble isoflavanoid structure"

    # Check for 3,4-dihydroderivative of coumestans
    coumestan_pattern = Chem.MolFromSmarts("[C@]12Oc3ccccc3[C@@]1(O)COc1ccccc21")
    if mol.HasSubstructMatch(coumestan_pattern):
        return False, "Appears to be a coumestan, not a pterocarpan"

    # Additional checks based on molecular properties
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 600:
        return False, "Molecular weight outside typical range for pterocarpans"

    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 2 or n_rotatable > 10:
        return False, "Number of rotatable bonds outside typical range for pterocarpans"

    return True, "Contains the pterocarpan scaffold with required fused rings, substitutions, and structural features"