"""
Classifies: CHEBI:16460 polyprenol phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol_phosphate(smiles: str):
    """
    Determines if a molecule is a polyprenol phosphate based on its SMILES string.
    A polyprenol phosphate is a polyprenol chain esterified with a phosphate group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        tuple: (bool, str). True if molecule is a polyprenol phosphate, False otherwise, along with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for phosphate group (allowing for deprotonation)
    phosphate_pattern = Chem.MolFromSmarts("[P](=[O])([O])([O,H])[O,H]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Define SMARTS for a repeating pattern of 3 isoprene units (head to tail)
    polyprenol_chain_pattern = Chem.MolFromSmarts("[C](C)([H])=[C][C][C][C](=[C][C]([H])[H])[C]([H])=[C][C][C][C](=[C][C]([H])[H])[C]([H])=[C][C]")
    chain_matches = mol.GetSubstructMatches(polyprenol_chain_pattern)
    
    if len(chain_matches) < 1: #minimum 1 chain of three isoprenes.
        return False, "Not enough repeating isoprene units in chain"


    # Check for connectivity of phosphate to the chain (at the end)
    connectivity_pattern = Chem.MolFromSmarts("[C]([H])=[C][C][C][C](=[C][C]([H])[H])[C]([H])=[C][C][C][C](=[C][C]([H])[H])[C]([H])=[C][C]O[P]")
    connectivity_matches = mol.GetSubstructMatches(connectivity_pattern)

    if len(connectivity_matches) < 1:
        return False, "Phosphate group is not connected to chain via oxygen at the end"
    
    #Check for saturated bonds in isoprenoid units
    saturated_pattern = Chem.MolFromSmarts("[C]([H])([H])[C][C][C][C]([H])([H])") #isoprenoid with saturated bonds.
    saturated_matches = mol.GetSubstructMatches(saturated_pattern)
    if len(saturated_matches) > 0:
          return False, "Chain contains saturated isoprenoid units"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
      return False, "Molecular weight too low for a polyprenol phosphate"

    # Check number of oxygen atoms - should not be a carbohydrate
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count > 5 :
        return False, "Too many oxygen atoms, likely a carbohydrate"
    

    return True, "Contains a polyprenol chain esterified with a phosphate group"