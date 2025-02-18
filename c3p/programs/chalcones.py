"""
Classifies: CHEBI:23086 chalcones
"""
"""
Classifies: CHEBI:23003 chalcone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    A chalcone is a ketone that is 1,3-diphenylpropenone (benzylideneacetophenone), ArCH=CH(=O)Ar, and its derivatives formed by substitution.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a chalcone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the core chalcone structure: Ar-CH=CH-C(=O)-Ar or Ar-CH2-CH2-C(=O)-Ar (dihydrochalcone)
    # More flexible patterns to account for substituted aromatic rings and heteroaromatic rings
    chalcone_pattern = Chem.MolFromSmarts("[a]1[a][a][a][a][a]1-[CH]=[CH]-[C](=O)-[a]2[a][a][a][a][a]2")
    dihydrochalcone_pattern = Chem.MolFromSmarts("[a]1[a][a][a][a][a]1-[CH2]-[CH2]-[C](=O)-[a]2[a][a][a][a][a]2")
    
    if not (mol.HasSubstructMatch(chalcone_pattern) or mol.HasSubstructMatch(dihydrochalcone_pattern)):
        return False, "No chalcone core structure found (Ar-CH=CH-C(=O)-Ar or Ar-CH2-CH2-C(=O)-Ar)"

    # Check for the presence of two aromatic rings (Ar)
    aromatic_rings = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetIsAromatic()]
    if len(aromatic_rings) < 2:
        return False, "Less than two aromatic rings found"

    # Check for the presence of the propenone chain (CH=CH-C=O or CH2-CH2-C=O)
    propenone_chain = Chem.MolFromSmarts("[CH]=[CH]-[C](=O)")
    dihydro_propenone_chain = Chem.MolFromSmarts("[CH2]-[CH2]-[C](=O)")
    
    if not (mol.HasSubstructMatch(propenone_chain) or mol.HasSubstructMatch(dihydro_propenone_chain)):
        return False, "No propenone chain (CH=CH-C=O or CH2-CH2-C=O) found"

    # Check for the presence of a ketone group (C=O)
    ketone_group = Chem.MolFromSmarts("[C](=O)")
    if not mol.HasSubstructMatch(ketone_group):
        return False, "No ketone group (C=O) found"

    # Additional checks to filter out non-chalcone molecules
    # Check molecular weight (chalcones typically have MW > 200)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for chalcone"

    # Check number of carbons (chalcones typically have > 10 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Too few carbons for chalcone"

    # Relaxed oxygen count check (chalcones can have more than 3 oxygens due to substitutions)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen found"

    # Check for additional functional groups that might disqualify the molecule
    # For example, chalcones should not have carboxylic acids, esters, or amides
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][#6]")
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Carboxylic acid group found, not a chalcone"
    if mol.HasSubstructMatch(ester_pattern):
        return False, "Ester group found, not a chalcone"
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Amide group found, not a chalcone"

    return True, "Contains the core chalcone structure (Ar-CH=CH-C(=O)-Ar or Ar-CH2-CH2-C(=O)-Ar) with possible substitutions on the aromatic rings"