"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: CHEBI:33716 dihydroflavonol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_dihydroflavonols(smiles: str):
    """
    Determines if a molecule is a dihydroflavonol based on its SMILES string.
    A dihydroflavonol is a hydroxyflavanone with a hydroxy group at position 3 of the heterocyclic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroflavonol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for the flavanone scaffold with a hydroxy group at position 3
    flavanone_pattern = Chem.MolFromSmarts("[O;R]1[C@@H]([C@H](C(=O)C2=C1C=CC=C2)O)[c;R]3ccc(O)cc3")
    if not mol.HasSubstructMatch(flavanone_pattern):
        return False, "Molecule does not contain the dihydroflavonol scaffold"
    
    # Count the number of hydroxy groups
    hydroxyl_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and atom.GetHybridization() == Chem.rdchem.HybridizationType.SP3)
    if hydroxyl_groups < 3:
        return False, "Molecule does not have enough hydroxy groups"
    
    # Check for additional double bonds or rings (to exclude flavones and flavonols)
    num_double_bonds = mol.GetNumBonds(Chem.BondType.DOUBLE)
    num_rings = mol.GetRingInfo().NumRings()
    if num_double_bonds > 3 or num_rings > 3:
        return False, "Molecule has too many double bonds or rings"
    
    # Check molecular weight - dihydroflavonols typically 200-500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 500:
        return False, "Molecular weight outside typical range for dihydroflavonols"
    
    return True, "Molecule meets the criteria for a dihydroflavonol"