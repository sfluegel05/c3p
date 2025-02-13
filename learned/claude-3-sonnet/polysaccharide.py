"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: CHEBI:18064 polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a biomacromolecule consisting of large numbers of monosaccharide
    residues linked glycosidically. This term is commonly used only for those containing
    more than ten monosaccharide residues.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for common monosaccharide units
    glucose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O")
    fructose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O")
    galactose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O")
    mannose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O")
    xylose_pattern = Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@@H](O)[C@@H]1O")

    # Count the number of monosaccharide units
    monosaccharide_count = sum(
        1 for pattern in [glucose_pattern, fructose_pattern, galactose_pattern, mannose_pattern, xylose_pattern]
        if mol.HasSubstructMatch(pattern)
    )

    # Define a general SMARTS pattern for glycosidic linkages
    glycosidic_linkage_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O)[C@H](O)[C@@H](O[C@@H]2[C@H](O)[C@H](O)[C@@H](O)[C@@H]2O)O[C@@H]1CO")

    # Count the number of glycosidic linkages
    glycosidic_linkage_count = len(mol.GetSubstructMatches(glycosidic_linkage_pattern))

    # Check if the molecule has a sufficient number of monosaccharide units and glycosidic linkages
    if monosaccharide_count >= 10 and glycosidic_linkage_count >= 10:
        # Calculate molecular weight
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)

        # Check if the molecular weight falls within a reasonable range for polysaccharides (500 - 1,000,000 Da)
        if 500 < mol_wt < 1000000:
            return True, "Contains more than 10 monosaccharide units linked by glycosidic bonds, and molecular weight is within the expected range for polysaccharides"
        else:
            return False, f"Molecular weight ({mol_wt:.2f} Da) is outside the expected range for polysaccharides"
    else:
        return False, f"Insufficient number of monosaccharide units ({monosaccharide_count}) or glycosidic linkages ({glycosidic_linkage_count})"