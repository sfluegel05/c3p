"""
Classifies: CHEBI:16412 lipopolysaccharide
"""
"""
Classifies: lipopolysaccharide
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_lipopolysaccharide(smiles: str):
    """
    Determines if a molecule is a lipopolysaccharide based on its SMILES string.
    A lipopolysaccharide consists of a lipid A moiety (disaccharide of glucosamine units with fatty acids),
    a core polysaccharide containing heptose and KDO units, and an O-antigen polysaccharide chain.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a lipopolysaccharide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glucosamine disaccharide (lipid A backbone)
    glucosamine_disaccharide = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@@H](N)[C@@H](O)[C@H]1O[C@H]2[C@@H](O)[C@@H](N)[C@@H](O)[C@H]2O")
    if not mol.HasSubstructMatch(glucosamine_disaccharide):
        return False, "No glucosamine disaccharide (lipid A backbone) found"

    # Look for 3-hydroxytetradecanoic acid chains attached to glucosamine units
    fatty_acid = Chem.MolFromSmarts("C(=O)CCCCCCCCCCCCC[OH]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid)
    if len(fatty_acid_matches) < 2:
        return False, f"Found {len(fatty_acid_matches)} fatty acid chains, need at least 2"

    # Look for heptose units (7-carbon sugar rings)
    heptose = Chem.MolFromSmarts("OC1[C@H](O)[C@@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H]1O")
    heptose_matches = mol.GetSubstructMatches(heptose)
    if len(heptose_matches) < 2:
        return False, f"Found {len(heptose_matches)} heptose units, need at least 2"

    # Look for KDO unit (8-carbon sugar acid)
    kdo = Chem.MolFromSmarts("O=C(C(O)C(O)C(O)C(O)C(O)C=O)CO")
    if not mol.HasSubstructMatch(kdo):
        return False, "No KDO unit found"

    # Check for polysaccharide chain (multiple glycosidic linkages)
    glycosidic_linkage = Chem.MolFromSmarts("[C@H]1(O[C@@H]2[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]2O)[O][C@@H](O)[C@H](O)[C@H](O)[C@H]1O")
    polysaccharide_matches = mol.GetSubstructMatches(glycosidic_linkage)
    if len(polysaccharide_matches) < 3:
        return False, f"Found {len(polysaccharide_matches)} glycosidic linkages, indicating insufficient polysaccharide chain"

    # Check for high molecular weight (LPS molecules are typically very large)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 1000:
        return False, "Molecular weight too low for lipopolysaccharide"

    # If all checks pass, classify as lipopolysaccharide
    return True, "Contains lipid A backbone with fatty acid chains, core polysaccharide with heptose and KDO units, and polysaccharide chains"

__metadata__ = {'chemical_class': {'name': 'lipopolysaccharide',
                                    'definition': 'Liposaccharide natural compounds consisting of a trisaccharide repeating unit (two heptose units and octulosonic acid) with oligosaccharide side chains and 3-hydroxytetradecanoic acid units (they are a major constituent of the cell walls of Gram-negative bacteria).'}}