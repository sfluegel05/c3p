"""
Classifies: CHEBI:39437 tocol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_tocol(smiles: str):
    """
    Determines if a molecule is a tocol based on its SMILES string.
    A tocol is a chromanol with a chroman-6-ol skeleton substituted at position 2 by a
    saturated or triply-unsaturated hydrocarbon chain consisting of three isoprenoid units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the chroman-6-ol core.
    chromanol_pattern = Chem.MolFromSmarts('[OX2]1[c]2[c]([c]([#6])[c]3[O][C]([#6])[C]3[c]2)[c]([#6])[c]1')
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "Chroman-6-ol core not found"
    
    # SMARTS pattern for a 15 carbon chain with single or double bonds (isoprenoid units) attached to position 2.
    isoprenoid_chain_pattern = Chem.MolFromSmarts('[#6][CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]')

    #Check that the position 2 carbon is linked to the isoprenoid chain, by verifying that there is a match to the isoprenoid chain in the mol.
    isoprenoid_matches = mol.GetSubstructMatches(isoprenoid_chain_pattern)
    if len(isoprenoid_matches) == 0:
       return False, "No isoprenoid chain found"

    # Check the number of carbon atoms in the isoprenoid chain (expecting 15)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    core_carbon_count = 9 # based on the core SMARTS, excluding the carbon attached to position 2
    isoprenoid_carbon_count = carbon_count - core_carbon_count
    if isoprenoid_carbon_count != 15:
       return False, f"Isoprenoid chain does not have 15 carbons. Found {isoprenoid_carbon_count} carbon atoms in the chain."

    #Check the molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300: # Tocols are generally above 300 Da
        return False, "Molecular weight too low for tocol"

    return True, "Contains chroman-6-ol core with a 15-carbon isoprenoid chain at position 2."