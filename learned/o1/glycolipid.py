"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: CHEBI:33563 glycolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    A glycolipid is a lipid molecule with a carbohydrate attached by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for a carbohydrate (sugar moiety)
    sugar_pattern = Chem.MolFromSmarts("C1(CO)OC(O)C(O)C(O)C1O")  # Simplified glucose ring
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No carbohydrate moiety detected"

    # Define SMARTS pattern for a lipid chain (long aliphatic chain)
    lipid_chain_pattern = Chem.MolFromSmarts("CCCCCCCC")  # Represents a hydrocarbon chain
    if not mol.HasSubstructMatch(lipid_chain_pattern):
        return False, "No lipid chain detected"

    # Check for glycosidic linkage between sugar and lipid
    # A glycosidic bond is an ether linkage connecting the sugar ring to another group
    glycosidic_bond_found = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        # Look for an oxygen atom connecting two carbons (C-O-C)
        if bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
            if (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6) or \
               (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8):
                # Check if one carbon is in sugar and the other is in lipid chain
                sugar_match = sugar_pattern.HasSubstructMatch(Chem.MolFromSmiles(atom1.GetSmarts())) or \
                              sugar_pattern.HasSubstructMatch(Chem.MolFromSmiles(atom2.GetSmarts()))
                lipid_match = lipid_chain_pattern.HasSubstructMatch(Chem.MolFromSmiles(atom1.GetSmarts())) or \
                              lipid_chain_pattern.HasSubstructMatch(Chem.MolFromSmiles(atom2.GetSmarts()))
                if sugar_match and lipid_match:
                    glycosidic_bond_found = True
                    break

    if not glycosidic_bond_found:
        return False, "No glycosidic linkage detected between lipid and sugar"

    return True, "Molecule contains lipid and carbohydrate moieties connected by a glycosidic bond"