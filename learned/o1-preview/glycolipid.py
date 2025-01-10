"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is defined as any member of class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage to a carbohydrate part (usually a mono-, di- or tri-saccharide).
    Some bacterial glycolipids have the sugar part acylated by one or more fatty acids and the glycerol part may be absent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for monosaccharide (simplified pattern matching a ring with oxygen and hydroxyls)
    monosaccharide_smarts = "[C;R]1([O;R])[C;R][C;R][C;R][C;R][O;R]1"
    monosaccharide = Chem.MolFromSmarts(monosaccharide_smarts)

    # Define SMARTS pattern for glycerol backbone with acyl groups at positions 1 and 2
    glycerol_diacyl_smarts = "OCC(OC(=O)[#6])[C;!R]OC(=O)[#6]"
    glycerol_diacyl = Chem.MolFromSmarts(glycerol_diacyl_smarts)

    # Define SMARTS pattern for glycosidic linkage at position 3 of glycerol
    glycosidic_linkage_smarts = "OCC(O[C;R])[C;!R]"
    glycosidic_linkage = Chem.MolFromSmarts(glycosidic_linkage_smarts)

    # Define SMARTS pattern for acylated sugar (bacterial glycolipids)
    acylated_sugar_smarts = "[C;R]([O;R])[O;R]C(=O)[#6]"
    acylated_sugar = Chem.MolFromSmarts(acylated_sugar_smarts)

    # Check for monosaccharide unit
    if mol.HasSubstructMatch(monosaccharide):
        # Check for glycerol backbone with diacyl groups
        if mol.HasSubstructMatch(glycerol_diacyl):
            # Check for glycosidic linkage between glycerol and sugar
            if mol.HasSubstructMatch(glycosidic_linkage):
                return True, "Contains 1,2-di-O-acylglycerol linked to carbohydrate via glycosidic bond"
            else:
                return False, "No glycosidic linkage between glycerol and carbohydrate found"
        else:
            # Check for acylated sugar (bacterial glycolipid)
            if mol.HasSubstructMatch(acylated_sugar):
                return True, "Contains acylated sugar moiety (bacterial glycolipid without glycerol)"
            else:
                return False, "No glycerol diacyl backbone or acylated sugar found"
    else:
        return False, "No carbohydrate moiety detected"

    return False, "Does not match glycolipid structural criteria"