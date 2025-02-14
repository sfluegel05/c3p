"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: CHEBI:38308 2,5-diketopiperazines
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_5_diketopiperazine(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is a piperazinone with a piperazine-2,5-dione skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for piperazine ring with two carbonyls at positions 2 and 5
    piperazine_pattern = Chem.MolFromSmarts("C1NC(=O)CN(C1=O)C")
    if not mol.HasSubstructMatch(piperazine_pattern):
        return False, "No piperazine-2,5-dione core found"

    # Check for reasonable number of rings (diketopiperazines usually have 1-3 rings)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 1 or n_rings > 3:
        return False, f"Unusual number of rings ({n_rings}) for a diketopiperazine"

    # Check for common functional groups like esters, ethers, amides
    allowed_fgroups = ["C(=O)O", "C-O-C", "C(=O)N"]
    fgroup_patterns = [Chem.MolFromSmarts(p) for p in allowed_fgroups]
    has_allowed_fgroup = any(mol.HasSubstructMatch(p) for p in fgroup_patterns)
    if not has_allowed_fgroup:
        return False, "No common functional groups found for diketopiperazines"

    # Check for reasonable molecular weight range (typically 200-800 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 800:
        return False, f"Molecular weight ({mol_wt:.2f} Da) outside typical range for diketopiperazines"

    # Enumerate stereoisomers and check core structure on each
    isomers = list(Chem.EnumerateStereoisomers(mol))
    for isomer in isomers:
        if isomer.HasSubstructMatch(piperazine_pattern):
            return True, "Contains piperazine-2,5-dione core with allowed substituents"

    return False, "No stereoisomer matched the piperazine-2,5-dione core"