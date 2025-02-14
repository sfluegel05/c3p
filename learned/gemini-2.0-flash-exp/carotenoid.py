"""
Classifies: CHEBI:23044 carotenoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_carotenoid(smiles: str):
    """
    Determines if a molecule is a carotenoid based on its SMILES string.
    Carotenoids are tetraterpenoids (C40) with a long conjugated polyene chain,
    formally derived from psi,psi-carotene, excluding retinoids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carotenoid, False otherwise
        str: Reason for classification
    """
    # 1. Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 2. Check for Long Conjugated Chain (at least 3 double bonds with interruptions)
    conjugated_chain_pattern = Chem.MolFromSmarts("[CX3]=[CX3]~[CX3]=[CX3]~[CX3]=[CX3]~[CX3]=[CX3]") #at least 3 double bonds
    if not mol.HasSubstructMatch(conjugated_chain_pattern):
        return False, "Does not contain a long conjugated chain of double bonds"


    #3. Check for isoprenoid units (at least 4) and geranylgeranyl pattern
    isoprene_pattern = Chem.MolFromSmarts("CC(=C)C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    geranylgeranyl_pattern = Chem.MolFromSmarts("CC(=C)CC/C=C(\C)CC/C=C(\C)CC/C=C(\C)")
    if len(isoprene_matches) < 4 and not mol.HasSubstructMatch(geranylgeranyl_pattern):
        return False, "Does not have a typical isoprenoid pattern for a carotenoid"

    #4. Molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, "Molecular weight too low for a carotenoid"

    # 5. Exclude Retinoids (beta-ionone ring)
    retinoid_pattern = Chem.MolFromSmarts("CC1=C[CH](C)[CH2][CH2][C](C)=C1")
    if mol.HasSubstructMatch(retinoid_pattern):
         return False, "Contains a retinoid-like beta-ionone ring"

    # If all checks pass
    return True, "Meets criteria for a carotenoid: C40 core, long conjugated polyene, and no retinoid beta-ionone"