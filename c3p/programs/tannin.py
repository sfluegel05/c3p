"""
Classifies: CHEBI:26848 tannin
"""
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

def is_tannin(smiles):
    """
    Determines if a molecule is a tannin based on its SMILES string.
    Tannins are polyphenolic compounds with a diverse range of structures, often containing galloyl, phenolic, catechol, pyrogallol, glucose, and ester subunits.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tannin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular properties
    mw = rdMolDescriptors.CalcExactMolWt(mol)
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    n_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    n_rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    
    # Count hydrogen bond donors and acceptors
    n_hydrogen_donors = sum(1 for atom in mol.GetAtoms() if atom.GetTotalNumHs() > 0)
    n_hydrogen_acceptors = sum(1 for atom in mol.GetAtoms() if atom.GetTotalNumHs() == 0 and atom.GetAtomicNum() in [7, 8, 16])
    
    # Look for common tannin substructures
    galloyl_pattern = Chem.MolFromSmarts("*c1c(O)c(O)c(O)c(C(=O)O)c1*")
    phenol_pattern = Chem.MolFromSmarts("*c1ccc(O)cc1*")
    catechol_pattern = Chem.MolFromSmarts("*c1c(O)c(O)ccc1*")
    pyrogallol_pattern = Chem.MolFromSmarts("*c1c(O)c(O)c(O)cc1*")
    glucose_pattern = Chem.MolFromSmarts("*OC1OC(CO)C(O)C(O)C1O*")
    ester_pattern = Chem.MolFromSmarts("*C(=O)O*")
    
    has_galloyl = mol.HasSubstructMatch(galloyl_pattern)
    has_phenol = mol.HasSubstructMatch(phenol_pattern)
    has_catechol = mol.HasSubstructMatch(catechol_pattern)
    has_pyrogallol = mol.HasSubstructMatch(pyrogallol_pattern)
    has_glucose = mol.HasSubstructMatch(glucose_pattern)
    has_ester = mol.HasSubstructMatch(ester_pattern)
    
    # Set up tannin classification criteria
    is_polyphenolic = n_hydrogen_donors >= 4 and n_hydrogen_acceptors >= 8
    has_tannin_substructures = has_galloyl or has_phenol or has_catechol or has_pyrogallol or has_glucose or has_ester
    complex_structure = n_rings >= 4 and n_aromatic_rings >= 2 and n_rotatable_bonds >= 6
    molecular_weight_range = 300 < mw < 3000
    
    if is_polyphenolic and has_tannin_substructures and complex_structure and molecular_weight_range:
        return True, "Molecule contains polyphenolic substructures and meets tannin criteria"
    else:
        return False, "Molecule does not meet tannin criteria"