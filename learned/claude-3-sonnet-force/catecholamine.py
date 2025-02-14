"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: CHEBI:18148 catecholamine
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    A catecholamine is a 4-(2-aminoethyl)pyrocatechol (4-(2-aminoethyl)benzene-1,2-diol)
    or a derivative formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pyrocatechol core with alkyl side chain containing amine
    catecholamine_pattern = Chem.MolFromSmarts("c1c(O)cc(O)cc1CCCN")
    catecholamine_matches = mol.GetSubstructMatches(catecholamine_pattern)
    if not catecholamine_matches:
        return False, "No catecholamine core found"

    # Check for additional hydroxyl groups or methyl substituents
    num_hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetDegree() == 1)
    num_methyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() == 4 and atom.GetTotalNumHs() == 3)

    # Check for quinone or oxidized form
    quinone_pattern = Chem.MolFromSmarts("c1c(=O)c(=O)cc(CCCN)c1")
    quinone_matches = mol.GetSubstructMatches(quinone_pattern)

    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 300:
        return False, "Molecular weight outside typical range for catecholamines"

    # Check for specific ring systems and atom counts
    ring_info = mol.GetRingInfo()
    if not any(len(ring) == 6 for ring in ring_info.AtomRings()):
        return False, "No benzene ring found"

    num_nitrogen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if num_nitrogen_atoms != 1:
        return False, "Must have exactly one nitrogen atom"

    num_oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if num_oxygen_atoms < 2 or num_oxygen_atoms > 4:
        return False, "Number of oxygen atoms outside typical range for catecholamines"

    # Classify based on structural features
    if quinone_matches:
        return True, "Oxidized catecholamine derivative (quinone form)"
    elif num_hydroxyl_groups > 2 or num_methyl_groups > 0:
        return True, "Substituted catecholamine derivative"
    else:
        return True, "Unsubstituted catecholamine"