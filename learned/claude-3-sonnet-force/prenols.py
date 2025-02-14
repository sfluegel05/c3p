"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: CHEBI:38732 prenols
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES string.
    Prenols are alcohols with the general formula H-[CH2C(Me)=CHCH2]nOH,
    where the carbon skeleton is composed of one or more isoprene units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX1H]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "No hydroxyl group found"
    
    # Check for isoprene units
    isoprene_pattern = Chem.MolFromSmarts("[CH2]=[C]([CH3])[CH]=[CH2]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) == 0:
        return False, "No isoprene units found"
    
    # Check for linear carbon skeleton
    carbon_skeleton = Chem.DeleteRedundantBonds(mol, Chem.RemoveHsFromMol(mol))
    sssr = Chem.GetSymmSSSR(carbon_skeleton)
    if len(sssr) > 1:
        return False, "Carbon skeleton not linear"
    
    # Check for terminal hydroxyl group
    hydroxyl_atom = next((atom for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and atom.GetTotalNumHs() == 1), None)
    if hydroxyl_atom is None:
        return False, "No terminal hydroxyl group found"
    
    # Check for isoprene units in linear skeleton
    isoprene_atoms = set([atom.GetIdx() for match in isoprene_matches for atom in match])
    skeleton_atoms = set([atom.GetIdx() for atom in carbon_skeleton.GetAtoms()])
    if not isoprene_atoms.issubset(skeleton_atoms):
        return False, "Isoprene units not part of linear carbon skeleton"
    
    # Check molecular weight and number of atoms
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    n_atoms = mol.GetNumAtoms()
    if mol_wt < 200 or n_atoms < 10:
        return False, "Molecule too small to be a prenol"
    
    return True, "Molecule contains a linear carbon skeleton with isoprene units and a terminal hydroxyl group"